import arcpy
import xml.etree.ElementTree as ET
import re

arcpy.ClearWorkspaceCache_management()


# ===============================================================
# QA TOOL: Cartographic Crossing Mask Detection
# ---------------------------------------------------------------
# Purpose:
#   Detects locations where a road or railway passes beneath
#   each other at a different vertical level (STUFE), and verifies
#   that a cartographic mask polygon exists at each such crossing.
#
# Pipeline overview (per process config):
#   1. Select input roads + crossing roads/rails from the GDB
#   2. Run a geometry (centre-line) intersect to find candidate
#      crossing pairs — this avoids generating masks for all features
#   3. Apply the STUFE rule: keep only pairs where the input road
#      is at a strictly lower level than the crossing feature
#   4. Generate solid Feature Outline Masks only for the surviving
#      features (dash effects stripped temporarily from the CIM)
#   5. Intersect the mask polygons to find visual overlap areas
#   6. For each overlap polygon, verify it represents a true crossing
#      using three conditions (see is_true_crossing)
#   7. Check whether existing masks in the GDB cover each confirmed
#      crossing, and flag missing or misplaced ones
#
# Inputs (from XML config per process):
#   - Input dataset + WHERE clause (the "lower" road/rail)
#   - Crossing dataset(s) + WHERE clause (the "upper" road/rail)
#   - AboveExpression: SQL WHERE clause expressing STUFE_input < STUFE_crossing
#   - MaskDataset: FC containing existing cartographic mask polygons
#   - MaskAssociation: table linking road features (LINEID) to masks (MASKID)
#
# Outputs:
#   - QA_ValidCrossings_<i>: per-config point FC, one row per crossing pair
#   - QA_ValidCrossings: merged result across all configs
#   Fields:
#     CROSS_SOURCE   — which crossing dataset produced this row
#     CROSS_TYPE     — "CROSSING" or "ENDPOINT_TOUCH" or "NO_POLYGON"
#     MASK_ERROR     — 0 = OK, 1 = error
#     MASK_ERROR_TYP — "OK", "NO_MASK", "MASK_MISPLACED", or "N/A"
# ===============================================================

# -------------------------------------------------------
# HELPERS
# -------------------------------------------------------


def parse_configuration(tree_xml):
    """
    Parse the CartographicProcess XML config file and extract all
    'CreateCrossingMasks' process definitions.

    Each process defines one crossing scenario.
    Returns a list of raw parameter dicts, one per matching process.
    """
    root = tree_xml.getroot()
    results = []
    for process in root.findall(".//Process"):
        type_ref = process.find("TypeReference")
        if type_ref is not None and type_ref.get("name") == "CreateCrossingMasks":
            name = process.get("name")

            def get_param(name):
                return [
                    p.get("value")
                    for p in process.findall(f".//Parameter[@name='{name}']")
                ]

            results.append(
                {
                    "Process": name,
                    "InputDataset": (
                        get_param("InputDataset")[0]
                        if get_param("InputDataset")
                        else None
                    ),
                    "CrossingDatasets": get_param("CrossingDatasets"),
                    "AboveExpression": (
                        get_param("AboveExpression")[0]
                        if get_param("AboveExpression")
                        else None
                    ),
                    "MaskDataset": (
                        get_param("MaskDataset")[0]
                        if get_param("MaskDataset")
                        else None
                    ),
                    "MaskAssociation": (
                        get_param("MaskAssociation")[0]
                        if get_param("MaskAssociation")
                        else None
                    ),
                }
            )
    return results


def build_config(raw_config, crossing_suffix="_1"):
    """
    Transform a raw parameter dict (from parse_configuration) into a
    structured config dict ready for the processing loop.

    Args:
        raw_config: dict returned by parse_configuration
        crossing_suffix: suffix applied to crossing-side fields after
                         Intersect (default "_1")

    Returns:
        Structured config dict with keys:
          input_ds, input_where, crossing_configs,
          AboveExpression, MaskDataset, MaskAssociation
    """

    def split_dataset_and_where(value):
        # Config values are semicolon-separated: "DATASET_NAME;WHERE_CLAUSE"
        if not value:
            return None, None
        parts = [p.strip() for p in value.split(";") if p.strip()]
        return parts[0], parts[1] if len(parts) > 1 else None

    def transform_expression(expr):
        # Rewrite "input.STUFE < crossing.STUFE" -> "STUFE < STUFE_1"
        # so it can be used directly as an ArcGIS SQL WHERE clause
        if not expr:
            return None
        expr = re.sub(r"input\.(\w+)", r"\1", expr)
        expr = re.sub(r"crossing\.(\w+)", rf"\1{crossing_suffix}", expr)
        return expr

    input_ds, input_where = split_dataset_and_where(raw_config.get("InputDataset"))
    crossing_configs = [
        {"ds": ds, "where": where}
        for ds, where in (
            split_dataset_and_where(c) for c in raw_config.get("CrossingDatasets", [])
        )
        if ds
    ]
    return {
        "Process": raw_config.get("Process"),
        "input_ds": input_ds,
        "input_where": input_where,
        "crossing_configs": crossing_configs,
        "AboveExpression": transform_expression(raw_config.get("AboveExpression")),
        "MaskDataset": raw_config.get("MaskDataset"),
        "MaskAssociation": raw_config.get("MaskAssociation"),
    }


def pt_to_ground_meters(width_pt, scale):
    return width_pt * (1 / 72) * 0.0254 * scale


def make_solid_cim(lyr):
    """
    Return a modified CIM (Cartographic Information Model) definition for
    the given layer with all dash/gap effects removed from every symbol class.

    This is necessary because FeatureOutlineMasks generates the polygon
    footprint of the rendered symbol. Dashed lines produce fragmented mask
    polygons with gaps, making it impossible to:
      - reliably count boundary intersection points, or
      - detect whether erasing an overlap polygon splits the mask in two.
    Generating a solid (continuous) mask avoids both problems.

    The original CIM is NOT modified — this function returns a new CIM
    object.

    Three types of dash effects are removed:
      1. CIMGeometricEffectDashes on the top-level symbol
      2. CIMGeometricEffectDashes on individual symbol layers
      3. dashTemplate properties on CIMSolidStroke layers
    """
    cim = lyr.getDefinition("V2")

    for group in cim.renderer.groups:
        for cls in group.classes:
            sym = cls.symbol.symbol

            # Top-level symbol effects (dashes applied above individual layers)
            if hasattr(sym, "effects") and sym.effects:
                sym.effects = [
                    e
                    for e in sym.effects
                    if e.__class__.__name__ != "CIMGeometricEffectDashes"
                ]

            # Per-layer effects and properties
            for sl in sym.symbolLayers:
                if hasattr(sl, "effects") and sl.effects:
                    sl.effects = [
                        e
                        for e in sl.effects
                        if e.__class__.__name__ != "CIMGeometricEffectDashes"
                    ]

                if hasattr(sl, "dashTemplate") and sl.dashTemplate:
                    sl.dashTemplate = None

                if sl.__class__.__name__ == "CIMVectorMarker":
                    if hasattr(sl, "enable"):
                        sl.enable = False

    return cim


def apply_solid_cim(lyr):
    """
    Apply the solid (dash-free) CIM to the given layer and return the
    original CIM so it can be restored afterwards.
    """
    original_cim = lyr.getDefinition("V2")
    lyr.setDefinition(make_solid_cim(lyr))
    return original_cim


def restore_cim(lyr, original_cim):
    lyr.setDefinition(original_cim)


def intermediate(name, i=None, debug=False):
    """
    Return a path for an intermediate FC.
    In debug mode: saves to GDB with a DEBUG_ prefix.
    In normal mode: saves to in_memory.
    """
    full_name = f"{name}_{i}" if i is not None else name
    if debug:
        return f"{gdb}/DEBUG_{full_name}"
    return f"in_memory/{full_name}"


def save_debug(fc, name, i, gdb):
    """
    Copy an in_memory FC to GDB for inspection (called only in debug mode).
    """
    out = f"{gdb}/DEBUG_{name}_{i}"
    if arcpy.Exists(out):
        arcpy.management.Delete(out)
    arcpy.management.CopyFeatures(fc, out)
    return out


def is_true_crossing(mask_geom, in_line_geo, cr_line_geo, overlap_polygon):
    """
    Determine whether the overlap polygon between two mask polygons
    represents a genuine visual crossing (not an endpoint touch).

    A true crossing occurs when one road's symbol physically passes
    through the other's symbol from side to side. An endpoint touch
    (T-junction or end-to-end connection) occurs when a road ends
    at or near the other road — these look like overlaps on the map
    but do not require a cartographic mask.

    For each part of the overlap polygon, two conditions must hold:

    CONDITION 1 — Erase splits the input mask:
        Erasing this overlap polygon part from the input mask polygon
        must increase the polygon's part count. If the mask splits into
        two separate pieces, the crossing road passed through it from
        one side to the other — a true crossing. If the mask stays in
        one piece, the overlap only nibbled an edge or end cap — not
        a true crossing.

    CONDITION 2 — No shared endpoint inside the overlap:
        If a point where the input and crossing lines share an endpoint
        (i.e. one road starts or ends where the other starts or ends)
        falls inside this overlap polygon part, the apparent "split" is
        caused by a road terminating at that location, not by a true
        crossing. Such parts are rejected.

    Args:
        mask_geom:       Polygon — the input road's solid mask polygon
        in_line_geo:     Polyline — the input road's centre-line geometry
        cr_line_geo:     Polyline — the crossing road's centre-line geometry
        overlap_polygon: Polygon — the intersection of the two mask polygons
                         (may be multipart if the same pair crosses multiple
                         times; each part is evaluated independently)

    Returns:
        (is_crossing: bool, valid_part_indices: list[int])
        valid_part_indices contains the part indices of overlap_polygon
        that passed both conditions; these are the confirmed crossing
        locations used for the mask QA check.
    """
    ok_parts = []
    try:
        sr_crossing = mask_geom.spatialReference
        in_endpoints = arcpy.Multipoint(
            arcpy.Array([in_line_geo.firstPoint, in_line_geo.lastPoint]), sr_crossing
        )
        cr_endpoints = arcpy.Multipoint(
            arcpy.Array([cr_line_geo.firstPoint, cr_line_geo.lastPoint]), sr_crossing
        )
        shared_endpoints = in_endpoints.intersect(cr_endpoints, 1)

        for p in range(overlap_polygon.partCount):
            part = arcpy.Polygon(overlap_polygon.getPart(p), sr_crossing)
            parts_before = mask_geom.partCount
            erased = mask_geom.difference(part)

            # CONDITION 1 — erasing this part must split the input mask
            if erased is None or erased.area == 0:
                continue
            if erased.partCount <= parts_before:
                continue

            # CONDITION 2 — no shared endpoint inside this part
            endpoint_overlap_in_inter = False
            if shared_endpoints:
                for pt in shared_endpoints:
                    pt_geom = arcpy.PointGeometry(pt, sr_crossing)
                    if not part.disjoint(pt_geom):
                        endpoint_overlap_in_inter = True
                        break

            if endpoint_overlap_in_inter:
                continue
            # Both conditions passed — this part is a confirmed crossing
            ok_parts.append(p)

    except Exception as e:
        print(f"  is_true_crossing ERROR: {e}")
    return len(ok_parts) > 0, ok_parts


def keep_parts(polygon, keep_indexes):
    """
    Return a new polygon containing only the parts at the given indices.
    Used after is_true_crossing to retain only the confirmed crossing
    polygons from the full overlap polygon, discarding endpoint-touch parts.
    The trimmed polygon is then used for the mask QA check so that a mask
    near an endpoint-touch location does not falsely satisfy the QA for a
    true crossing at a different location on the same road pair.
    """
    if not keep_indexes:
        return polygon
    parts = arcpy.Array()
    for idx in keep_indexes:
        parts.add(polygon.getPart(idx))
    return arcpy.Polygon(parts, polygon.spatialReference)


# -------------------------------------------------------
# CONFIG
# -------------------------------------------------------

# Set DEBUG = True to save all intermediate FCs to the GDB with
# a DEBUG_ prefix for inspection in ArcGIS Pro.
# Set DEBUG = False for production runs (uses in_memory only).
DEBUG = True

tree = ET.parse(r"J:\System\K2\20_Konfig\DKM25\DKM25_CartoProcess_Config.cp.xml")

gdb = r"D:\K2\WorkUnitData\WU1207_Thun\WU1207_Thun.gdb"
aprx_path = r"D:\K2\WorkUnitData\WU1207_Thun\WU1207_Thun.aprx"
aprx = arcpy.mp.ArcGISProject(aprx_path)
m = aprx.listMaps()[0]

arcpy.env.workspace = gdb
arcpy.env.overwriteOutput = True

configs = parse_configuration(tree)
config = build_config(
    configs[2]
)  # Change the index to run a different process definition.
input_ds = config["input_ds"]
input_where = config["input_where"]
crossing_configs = config["crossing_configs"]
above = config["AboveExpression"]
mask_fc = config["MaskDataset"]
mask_tb = config["MaskAssociation"]
sr = arcpy.Describe(input_ds).spatialReference
reference_scale = 25000
error_fc = "QA_ValidCrossings"
in_fid_field = f"FID_{input_ds}"

print(f"Process: {config['Process']}")
print(f"Debug mode: {DEBUG}")

# -------------------------------------------------------
# CLEANUP
# -------------------------------------------------------
# Remove any layers and GDB FCs left over from a previous run.
# Also clears in_memory workspace and compacts the GDB to release
# any lingering schema locks that would block FeatureOutlineMasks.

print("CLEANUP")
prefixes = ("DEBUG_", "QA_")
for lyr in list(m.listLayers()):
    try:
        if lyr.name.startswith(prefixes):
            m.removeLayer(lyr)
    except Exception as e:
        print(f"  WARNING layer remove failed: {lyr.name} -> {e}")
aprx.save()  # save after removing layers so locks are released

arcpy.env.workspace = gdb
datasets = []
for fc in arcpy.ListFeatureClasses():
    if fc.startswith(prefixes):
        datasets.append(fc)
for tb in arcpy.ListTables():
    if tb.startswith(prefixes):
        datasets.append(tb)
for ds in datasets:
    try:
        arcpy.management.Delete(ds)
    except Exception as e:
        print(f"  WARNING delete failed: {ds} -> {e}")

try:
    arcpy.env.workspace = "in_memory"
    mem_fcs = arcpy.ListFeatureClasses() or []
    mem_tables = arcpy.ListTables() or []
    for ds in mem_fcs + mem_tables:
        try:
            arcpy.management.Delete(ds)
        except Exception as e:
            print(f"  WARNING memory delete failed: {ds} -> {e}")
except Exception as e:
    print(f"  WARNING in_memory cleanup failed: {e}")

arcpy.env.workspace = gdb
arcpy.ClearWorkspaceCache_management()

try:
    arcpy.management.Compact(gdb)
except Exception as e:
    print(f"  WARNING compact failed: {e}")
print("  Done")

# -------------------------------------------------------
# 1. SELECT INPUTS + GENERATE INPUT MASKS
# -------------------------------------------------------
# The input mask is generated once and shared across all crossing configs
# in the loop below. This avoids regenerating the same input masks
# once per crossing dataset.

print("SELECT INPUTS + GENERATE INPUT MASKS")
input_aprx_lyr = m.listLayers(input_ds)[0]
arcpy.management.SelectLayerByAttribute(input_aprx_lyr, "NEW_SELECTION", input_where)
print(f"  Inputs: {int(arcpy.management.GetCount(input_aprx_lyr).getOutput(0))}")

# Temporarily make the symbol solid (no dashes) so FeatureOutlineMasks
# produces a continuous polygon. The original CIM is restored in the
# finally block regardless of whether the tool succeeds or fails.
original_input_cim = apply_solid_cim(input_aprx_lyr)
aprx.save()  # save is required for the CIM change to be visible to the tool

try:
    arcpy.cartography.FeatureOutlineMasks(
        input_aprx_lyr,
        f"{gdb}/QA_InputMasks",
        reference_scale,
        sr,
        "0 Meters",
        "EXACT",
        "ALL_FEATURES",
        "ALL",
    )
finally:
    restore_cim(input_aprx_lyr, original_input_cim)
    aprx.save()

print(
    f"  Input masks: {int(arcpy.management.GetCount(f'{gdb}/QA_InputMasks').getOutput(0))}"
)
arcpy.management.SelectLayerByAttribute(input_aprx_lyr, "CLEAR_SELECTION")

input_mask_geom = {}
with arcpy.da.SearchCursor(f"{gdb}/QA_InputMasks", [in_fid_field, "SHAPE@"]) as cur:
    for fid, geom in cur:
        input_mask_geom[fid] = geom

input_line_geom = {}
with arcpy.da.SearchCursor(
    input_ds, ["OBJECTID", "SHAPE@"], where_clause=input_where
) as cur:
    for oid, geom in cur:
        input_line_geom[oid] = geom

# -------------------------------------------------------
# 2. LOAD QA MASK TABLES
# -------------------------------------------------------
# Load the existing cartographic mask data once, shared across all configs.
#
# mask_lookup:      LINEID -> [MASKID, ...]
#   Maps each road feature (by its OBJECTID as LINEID) to the list of
#   existing mask polygon IDs associated with it.
#
# mask_geom_lookup: MASKID -> Polygon geometry
#   Maps each existing mask polygon OBJECTID to its geometry, used
#   to check whether a mask actually overlaps the crossing location.

print("LOAD QA MASK TABLES")
mask_lookup = {}
with arcpy.da.SearchCursor(mask_tb, ["LINEID", "MASKID"]) as cur:
    for lineid, maskid in cur:
        mask_lookup.setdefault(lineid, []).append(maskid)

mask_geom_lookup = {}
with arcpy.da.SearchCursor(mask_fc, ["OBJECTID", "SHAPE@"]) as cur:
    for oid, geom in cur:
        mask_geom_lookup[oid] = geom

print(f"  Lines with masks: {len(mask_lookup)}")
print(f"  Mask geometries:  {len(mask_geom_lookup)}")

# -------------------------------------------------------
# 3. PER-CONFIG PROCESSING LOOP
# -------------------------------------------------------
# Each crossing config defines one crossing dataset (e.g. all roads,
# or all railways). The loop processes them independently and collects
# one output FC per config, merged at the end.

result_fcs = []

for i, cfg in enumerate(crossing_configs):
    print(f"\nPROCESSING CONFIG {i}: {cfg['ds']}")

    cross_ds = cfg["ds"]
    cross_where = cfg["where"]

    # ---------------------------------------------------
    # 3.1 SELECT CROSSINGS + GENERATE CROSSING MASKS
    # ---------------------------------------------------
    # Same approach as the input masks: select features, strip dashes,
    # run FeatureOutlineMasks, restore CIM immediately after.

    cross_aprx_lyr = m.listLayers(cross_ds)[0]
    arcpy.management.SelectLayerByAttribute(
        cross_aprx_lyr, "NEW_SELECTION", cross_where
    )
    print(f"  Crossings: {int(arcpy.management.GetCount(cross_aprx_lyr).getOutput(0))}")

    cross_mask_fc = f"{gdb}/QA_CrossMasks_{i}"
    original_cross_cim = apply_solid_cim(cross_aprx_lyr)
    aprx.save()

    try:
        arcpy.cartography.FeatureOutlineMasks(
            cross_aprx_lyr,
            cross_mask_fc,
            reference_scale,
            sr,
            "0 Meters",
            "EXACT",
            "ALL_FEATURES",
            "ALL",
        )
    finally:
        restore_cim(cross_aprx_lyr, original_cross_cim)
        aprx.save()

    print(
        f"  Cross masks: {int(arcpy.management.GetCount(cross_mask_fc).getOutput(0))}"
    )
    arcpy.management.SelectLayerByAttribute(cross_aprx_lyr, "CLEAR_SELECTION")

    # FeatureOutlineMasks writes a FID_<dataset> field into the mask FC
    # linking each mask polygon to its source feature's OBJECTID.
    # Detect the field name dynamically since dataset names vary.
    cr_fid_field = next(
        (
            f.name
            for f in arcpy.ListFields(cross_mask_fc)
            if f.name.startswith("FID_DKM25")
        ),
        None,
    )
    if cr_fid_field is None:
        print(f"  WARNING: cr_fid field not found, skipping config {i}")
        continue

    cross_mask_geom = {}
    with arcpy.da.SearchCursor(cross_mask_fc, [cr_fid_field, "SHAPE@"]) as cur:
        for fid, geom in cur:
            cross_mask_geom[fid] = geom

    cross_line_geom = {}
    with arcpy.da.SearchCursor(
        cross_ds, ["OBJECTID", "SHAPE@"], where_clause=cross_where
    ) as cur:
        for oid, geom in cur:
            cross_line_geom[oid] = geom

    # ---------------------------------------------------
    # 3.2 INTERSECT MASKS
    # ---------------------------------------------------
    # Two separate intersects are run:
    #
    # POLYGON output: the area where the input and crossing mask polygons
    #   overlap. This is the "overlap polygon" used in is_true_crossing
    #   to test whether erasing it from the input mask splits it in two.
    #
    # POINT output: the boundary intersection points between the two mask
    #   polygon outlines. The number of these points is used as a fast
    #   pre-filter — a true crossing requires the crossing line to enter
    #   and exit both sides of the input mask, producing at least 4
    #   boundary contact points. Endpoint touches produce 1–3 points.

    inter_poly_path = intermediate("inter_poly", i, DEBUG)
    inter_pt_path = intermediate("inter_pt", i, DEBUG)

    inter_poly = arcpy.analysis.Intersect(
        [f"{gdb}/QA_InputMasks", cross_mask_fc],
        inter_poly_path,
        "ALL",
        output_type="INPUT",
    )
    inter_pt = arcpy.analysis.Intersect(
        [f"{gdb}/QA_InputMasks", cross_mask_fc],
        inter_pt_path,
        "ALL",
        output_type="POINT",
    )

    inter_cr_fid = next(
        (
            f.name
            for f in arcpy.ListFields(inter_pt)
            if f.name.startswith("FID_DKM25") and f.name != in_fid_field
        ),
        None,
    )
    if inter_cr_fid is None:
        print(f"  WARNING: inter_cr_fid not found, skipping config {i}")
        continue

    print(
        f"  Intersection points: {int(arcpy.management.GetCount(inter_pt).getOutput(0))}"
    )

    # ---------------------------------------------------
    # 3.3 APPLY ABOVE RULE + DISSOLVE
    # ---------------------------------------------------
    # Filter both the point and polygon intersect results to keep only
    # pairs where STUFE_input < STUFE_crossing (the "above" expression
    # from the XML config). This eliminates same-level intersections and
    # cases where the input road is actually above the crossing feature.
    #
    # Dissolve to one row per (input_fid, crossing_fid) pair so that
    # each unique road combination is evaluated exactly once regardless
    # of how many intersection points or polygon fragments it produced.
    # MULTI_PART keeps all fragments as a single multipart geometry.

    sel_pt_path = intermediate("sel_pt", i, DEBUG)
    sel_poly_path = intermediate("sel_poly", i, DEBUG)

    sel_pt = arcpy.analysis.Select(inter_pt, sel_pt_path, above)
    sel_poly = arcpy.analysis.Select(inter_poly, sel_poly_path, above)
    print(
        f"  After above filter: {int(arcpy.management.GetCount(sel_pt).getOutput(0))}"
    )

    dissolve_fields = [in_fid_field, inter_cr_fid]

    valid_pt_path = intermediate("valid_pt", i, DEBUG)
    valid_poly_path = intermediate("valid_poly", i, DEBUG)

    valid_pt = arcpy.management.Dissolve(
        sel_pt, valid_pt_path, dissolve_fields, multi_part="MULTI_PART"
    )
    valid_poly = arcpy.management.Dissolve(
        sel_poly, valid_poly_path, dissolve_fields, multi_part="MULTI_PART"
    )

    arcpy.management.AddField(valid_pt, "ORIG_OID", "LONG")
    with arcpy.da.UpdateCursor(valid_pt, ["OBJECTID", "ORIG_OID"]) as cur:
        for row in cur:
            row[1] = row[0]
            cur.updateRow(row)

    # ---------------------------------------------------
    # 3.4 BUILD LOOKUPS
    # ---------------------------------------------------
    # Build in-memory dictionaries for fast per-row lookup in the
    # detection loop. All geometry operations happen in Python memory
    # without further tool calls, which is significantly faster than
    # running spatial analysis tools per crossing pair.

    # overlap_poly_lookup: (in_fid, cr_fid) -> overlap polygon geometry
    # If the same pair produces multiple polygon fragments (e.g. two
    # separate crossing locations), union them into one multipart polygon
    # so all parts are evaluated together in is_true_crossing.

    overlap_poly_lookup = {}
    with arcpy.da.SearchCursor(
        valid_poly, [in_fid_field, inter_cr_fid, "SHAPE@"]
    ) as cur:
        for in_fid, cr_fid, poly_geom in cur:
            if cr_fid is not None and cr_fid != -1:
                key = (in_fid, cr_fid)
                overlap_poly_lookup[key] = (
                    overlap_poly_lookup[key].union(poly_geom)
                    if key in overlap_poly_lookup
                    else poly_geom
                )

    pt_count_lookup = {}
    with arcpy.da.SearchCursor(valid_pt, [in_fid_field, inter_cr_fid, "SHAPE@"]) as cur:
        for in_fid, cr_fid, geom in cur:
            if cr_fid is not None and cr_fid != -1 and geom is not None:
                pt_count_lookup[(in_fid, cr_fid)] = geom.partCount

    # ---------------------------------------------------
    # 3.5 TRUE CROSSING FILTER
    # ---------------------------------------------------
    # Classify each candidate crossing pair as a true crossing or an
    # endpoint touch using a three-stage filter:
    #
    # Stage 1 (no polygon): if no overlap polygon exists for the pair,
    #   the mask contact is too small to determine — treat as invalid.
    #
    # Stage 2 (point count): if fewer than 4 boundary contact points
    #   exist, the contact cannot be a true crossing. A true crossing
    #   requires the crossing line to enter and exit both sides of the
    #   input mask, which always produces >= 4 boundary points.
    #   Endpoint touches (T-junctions, end-to-end) produce 1–3 points.
    #
    # Stage 3 (is_true_crossing): for pairs that pass stages 1 and 2,
    #   apply the geometric erase-and-split test (see is_true_crossing).
    #   This handles complex cases where 4+ points exist but the contact
    #   is still not a true crossing (e.g. a road with a rounded end cap
    #   that produces 4 boundary points at the end of the mask).

    endpoint_oids = set()
    no_polygon_oids = set()
    true_cross_oids = set()

    with arcpy.da.SearchCursor(
        valid_pt, ["ORIG_OID", in_fid_field, inter_cr_fid]
    ) as cur:
        for orig_oid, in_fid, cross_fid in cur:
            if cross_fid is None or cross_fid == -1:
                continue

            key = (in_fid, cross_fid)
            overlap_poly = overlap_poly_lookup.get(key)

            if not overlap_poly:
                no_polygon_oids.add(orig_oid)
                continue

            if pt_count_lookup.get(key, 0) < 4:
                endpoint_oids.add(orig_oid)
                continue

            in_mask = input_mask_geom.get(in_fid)
            in_line = input_line_geom.get(in_fid)
            cr_line = cross_line_geom.get(cross_fid)

            crossing, valid_parts = (
                is_true_crossing(in_mask, in_line, cr_line, overlap_poly)
                if in_mask and in_line and cr_line
                else (False, [])
            )

            if crossing and valid_parts:
                overlap_poly_lookup[key] = keep_parts(overlap_poly, valid_parts)
                true_cross_oids.add(orig_oid)
            else:
                endpoint_oids.add(orig_oid)

    print(f"  No polygon:       {len(no_polygon_oids)}")
    print(f"  Endpoint touches: {len(endpoint_oids)}")
    print(f"  True crossings:   {len(true_cross_oids)}")

    # ---------------------------------------------------
    # 3.6 BUILD OUTPUT
    # ---------------------------------------------------
    # Write one point feature per (input road, crossing road) pair.
    # The point geometry is the dissolved multipoint from the mask
    # boundary intersection.
    #
    # Output fields:
    #   CROSS_SOURCE:   name of the crossing dataset for this config
    #   CROSS_TYPE:     classification result from stage 3.5
    #   MASK_ERROR:     0 = no error, 1 = mask missing or misplaced
    #   MASK_ERROR_TYP: "OK", "NO_MASK", "MASK_MISPLACED", or "N/A"
    #
    # Mask QA logic (for confirmed crossings):
    #   Look up all mask IDs associated with the input road via
    #   mask_lookup (LINEID -> [MASKID]). For each confirmed crossing
    #   polygon part, check whether any of those mask polygons overlaps
    #   it (not disjoint). ALL parts must be covered — if a road crosses
    #   the same road twice, both locations need a mask.

    out_fc_i = f"{gdb}/QA_ValidCrossings_{i}"

    needs_create = not arcpy.Exists(out_fc_i)
    if needs_create:
        arcpy.management.CopyFeatures(valid_pt, out_fc_i)
        for fname, ftype, flen in [
            ("CROSS_SOURCE", "TEXT", 30),
            ("CROSS_TYPE", "TEXT", 20),
            ("MASK_ERROR", "SHORT", 0),
            ("MASK_ERROR_TYP", "TEXT", 50),
        ]:
            arcpy.management.AddField(
                out_fc_i, fname, ftype, field_length=flen if flen else None
            )
    else:
        arcpy.management.DeleteRows(out_fc_i)
        arcpy.management.Append(valid_pt, out_fc_i, "NO_TEST")

    update_fields = [
        "SHAPE@",
        "ORIG_OID",
        in_fid_field,
        inter_cr_fid,
        "CROSS_SOURCE",
        "CROSS_TYPE",
        "MASK_ERROR",
        "MASK_ERROR_TYP",
    ]

    with arcpy.da.UpdateCursor(out_fc_i, update_fields) as cur:
        for row in cur:
            oid = row[1]
            line_fid = row[2]
            cross_fid = row[3]
            row[4] = cross_ds

            if oid in endpoint_oids or oid in no_polygon_oids:
                row[5] = "ENDPOINT_TOUCH" if oid in endpoint_oids else "NO_POLYGON"
                row[6] = 0
                row[7] = "N/A"
            elif line_fid not in mask_lookup:
                row[5] = "CROSSING"
                row[6] = 1
                row[7] = "NO_MASK"
            else:
                row[5] = "CROSSING"
                valid_poly = overlap_poly_lookup.get((line_fid, cross_fid))
                if valid_poly is None:
                    row[6] = 0
                    row[7] = "N/A"
                else:
                    all_covered = all(
                        any(
                            mask_geom_lookup.get(mid) is not None
                            and not mask_geom_lookup[mid].disjoint(
                                arcpy.Polygon(valid_poly.getPart(p), sr)
                            )
                            for mid in mask_lookup[line_fid]
                        )
                        for p in range(valid_poly.partCount)
                    )
                    row[6] = 0 if all_covered else 1
                    row[7] = "OK" if all_covered else "MASK_MISPLACED"

            cur.updateRow(row)

    result_fcs.append(out_fc_i)
    print(f"  Output: {int(arcpy.management.GetCount(out_fc_i).getOutput(0))} records")

# -------------------------------------------------------
# MERGE + SUMMARY
# -------------------------------------------------------
# Combine all per-config output FCs into a single QA_ValidCrossings FC.
# Uses DeleteRows + Append on subsequent runs to preserve the FC
# (and its symbology in ArcGIS Pro) rather than recreating it.

print("\nMERGE RESULTS")
out_fc = f"{gdb}/{error_fc}"
for lyr in m.listLayers(error_fc):
    m.removeLayer(lyr)

if arcpy.Exists(out_fc):
    arcpy.management.DeleteRows(out_fc)
    arcpy.management.Append(result_fcs, out_fc, "NO_TEST")
else:
    arcpy.management.Merge(result_fcs, out_fc)

print("\nSUMMARY")
counts = {"OK": 0, "NO_MASK": 0, "MASK_MISPLACED": 0, "N/A": 0}
with arcpy.da.SearchCursor(out_fc, ["MASK_ERROR_TYP"]) as cur:
    for (typ,) in cur:
        if typ in counts:
            counts[typ] += 1

total = sum(counts.values())
print(f"  Total records:          {total}")
print(f"  Endpoint touches (N/A): {counts['N/A']}")
print(f"  OK (mask present):      {counts['OK']}")
print(f"  No mask at all:         {counts['NO_MASK']}")
print(f"  Mask misplaced:         {counts['MASK_MISPLACED']}")
print("DONE")
