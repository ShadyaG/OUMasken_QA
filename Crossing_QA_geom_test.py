import arcpy
from arcpy.cim import CIMSolidStroke
import xml.etree.ElementTree as ET
import re

arcpy.ClearWorkspaceCache_management()

# -------------------------------------------------------
# HELPERS
# -------------------------------------------------------


def parse_configuration(tree_xml):
    root = tree_xml.getroot()
    results = []
    for process in root.findall(".//Process"):
        type_ref = process.find("TypeReference")
        if type_ref is not None and type_ref.get("name") == "CreateCrossingMasks":
            name = process.get("name")

            def get_param(n):
                return [
                    p.get("value")
                    for p in process.findall(f".//Parameter[@name='{n}']")
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
    def split_dataset_and_where(value):
        if not value:
            return None, None
        parts = [p.strip() for p in value.split(";") if p.strip()]
        return parts[0], parts[1] if len(parts) > 1 else None

    def transform_expression(expr):
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
    original_cim = lyr.getDefinition("V2")
    lyr.setDefinition(make_solid_cim(lyr))
    return original_cim


def restore_cim(lyr, original_cim):
    lyr.setDefinition(original_cim)


def intermediate(name, i=None, debug=False):
    full_name = f"{name}_{i}" if i is not None else name
    return f"{gdb}/DEBUG_{full_name}" if debug else f"in_memory/{full_name}"


def read_symbol_widths(layer_name, m, reference_scale):
    lyr = m.listLayers(layer_name)[0]
    cim = lyr.getDefinition("V2")
    widths = {}
    for group in cim.renderer.groups:
        for cls in group.classes:
            symbol_value = int(cls.values[0].fieldValues[0])
            max_width = 0
            for sl in cls.symbol.symbol.symbolLayers:
                if isinstance(sl, CIMSolidStroke):
                    max_width = max(max_width, sl.width)
            if max_width > 0:
                widths[symbol_value] = pt_to_ground_meters(max_width, reference_scale)
    return widths


def is_true_crossing(
    mask_geom, in_line_geom, cr_line_geom, overlap_poly, geom_inter_pts
):
    """
    For each part of the overlap polygon, a valid crossing requires:
      1. The part contains at least one geometry intersection point
      2. Erasing the part from the input mask increases part count
      3. No shared endpoint of input+crossing lines falls inside the part
    Returns (is_crossing, list_of_valid_part_indices).
    """
    valid_parts = []
    try:
        sr = mask_geom.spatialReference

        # Shared endpoints between the two lines
        in_endpoints = arcpy.Multipoint(
            arcpy.Array([in_line_geom.firstPoint, in_line_geom.lastPoint]), sr
        )
        cr_endpoints = arcpy.Multipoint(
            arcpy.Array([cr_line_geom.firstPoint, cr_line_geom.lastPoint]), sr
        )
        shared_endpoints = in_endpoints.intersect(cr_endpoints, 1)

        for p in range(overlap_poly.partCount):
            part = arcpy.Polygon(overlap_poly.getPart(p), sr)

            # CONDITION 1 — part must contain at least one geometry intersection point
            part_has_geom_inter = any(
                not part.disjoint(arcpy.PointGeometry(pt, sr)) for pt in geom_inter_pts
            )
            if not part_has_geom_inter:
                continue

            # CONDITION 2 — erasing this part must split the input mask
            parts_before = mask_geom.partCount
            erased = mask_geom.difference(part)
            if erased is None or erased.area == 0:
                continue
            if erased.partCount <= parts_before:
                continue

            # CONDITION 3 — no shared endpoint inside this part
            if shared_endpoints:
                endpoint_in_part = any(
                    not part.disjoint(arcpy.PointGeometry(pt, sr))
                    for pt in shared_endpoints
                )
                if endpoint_in_part:
                    continue

            valid_parts.append(p)

    except Exception as e:
        print(f"  is_true_crossing ERROR: {e}")

    return len(valid_parts) > 0, valid_parts


def keep_parts(polygon, keep_indexes):
    if not keep_indexes:
        return polygon
    parts = arcpy.Array()
    for idx in keep_indexes:
        parts.add(polygon.getPart(idx))
    return arcpy.Polygon(parts, polygon.spatialReference)


# -------------------------------------------------------
# CONFIG
# -------------------------------------------------------

DEBUG = True

tree = ET.parse(r"J:\System\K2\20_Konfig\DKM25\DKM25_CartoProcess_Config.cp.xml")

gdb = r"D:\K2\WorkUnitData\WU1207_Thun_new\WU1207_Thun.gdb"
aprx_path = r"D:\K2\WorkUnitData\WU1207_Thun_new\WU1207_Thun.aprx"
aprx = arcpy.mp.ArcGISProject(aprx_path)
m = aprx.listMaps()[0]

arcpy.env.workspace = gdb
arcpy.env.overwriteOutput = True

configs = parse_configuration(tree)
config = build_config(configs[2])
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

print("CLEANUP")
prefixes = ("DEBUG_", "QA_")
for lyr in list(m.listLayers()):
    try:
        if lyr.name.startswith(prefixes):
            m.removeLayer(lyr)
    except Exception as e:
        print(f"  WARNING layer remove failed: {lyr.name} -> {e}")
aprx.save()

arcpy.env.workspace = gdb
for fc in arcpy.ListFeatureClasses():
    if fc.startswith(prefixes):
        try:
            arcpy.management.Delete(fc)
        except Exception as e:
            print(f"  WARNING delete failed: {fc} -> {e}")

arcpy.env.workspace = "in_memory"
for ds in (arcpy.ListFeatureClasses() or []) + (arcpy.ListTables() or []):
    try:
        arcpy.management.Delete(ds)
    except:
        pass

arcpy.env.workspace = gdb
arcpy.ClearWorkspaceCache_management()
arcpy.management.Compact(gdb)
print("  Done")

# -------------------------------------------------------
# READ SYMBOL WIDTHS
# -------------------------------------------------------

print("READ CIM SYMBOLOGY")
input_symbol_widths = read_symbol_widths(input_ds, m, reference_scale)
cross_symbol_widths = {}
for cfg in crossing_configs:
    cross_symbol_widths.update(read_symbol_widths(cfg["ds"], m, reference_scale))
print(f"  Input symbols: {len(input_symbol_widths)}")
print(f"  Cross symbols: {len(cross_symbol_widths)}")

# -------------------------------------------------------
# LOAD QA MASK TABLES (once, shared across configs)
# -------------------------------------------------------

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
# PER-CONFIG PROCESSING LOOP
# -------------------------------------------------------

result_fcs = []

for i, cfg in enumerate(crossing_configs):
    print(f"\n{'='*60}")
    print(f"PROCESSING CONFIG {i}: {cfg['ds']}")
    print(f"{'='*60}")
    # -------------------------------------------------------
    # 1. SELECT INPUTs
    # -------------------------------------------------------

    print("\nSELECT INPUT")
    input_aprx_lyr = m.listLayers(input_ds)[0]
    arcpy.management.SelectLayerByAttribute(
        input_aprx_lyr, "NEW_SELECTION", input_where
    )
    print(
        f"  Input roads: "
        f"{int(arcpy.management.GetCount(input_aprx_lyr).getOutput(0))}"
    )

    # ---------------------------------------------------
    # 2. SELECT CROSSING FEATURES
    # ---------------------------------------------------
    cross_ds = cfg["ds"]
    cross_where = cfg["where"]

    print("  SELECT CROSSING FEATURES")
    cross_aprx_lyr = m.listLayers(cross_ds)[0]

    arcpy.management.SelectLayerByAttribute(
        cross_aprx_lyr, "NEW_SELECTION", cross_where
    )
    print(
        f"  Crossings: "
        f"{int(arcpy.management.GetCount(cross_aprx_lyr).getOutput(0))}"
    )

    # ---------------------------------------------------
    # 3. GEOMETRY INTERSECT — find pairs that actually cross
    # ---------------------------------------------------

    print("  GEOMETRY INTERSECT")

    geom_inter_pt = arcpy.analysis.Intersect(
        [input_aprx_lyr, cross_aprx_lyr],
        intermediate("geom_inter_pt", i, DEBUG),
        "ALL",
        output_type="POINT",
    )

    print(
        f"  Geometry intersection points (raw): "
        f"{int(arcpy.management.GetCount(geom_inter_pt).getOutput(0))}"
    )
    # Clear all
    arcpy.management.SelectLayerByAttribute(input_aprx_lyr, "CLEAR_SELECTION")
    arcpy.management.SelectLayerByAttribute(cross_aprx_lyr, "CLEAR_SELECTION")

    # Confirm field names produced by the intersect
    inter_cr_fid = next(
        (
            f.name
            for f in arcpy.ListFields(geom_inter_pt)
            if f.name.startswith("FID_DKM25") and f.name != in_fid_field
        ),
        None,
    )
    if inter_cr_fid is None:
        print(f"  WARNING: inter_cr_fid not found, skipping config {i}")
        continue

    # ---------------------------------------------------
    # 4. STUFE RULE on geometry intersect
    # ---------------------------------------------------
    print("  APPLY STUFE RULE")

    sel_geom_pt = arcpy.analysis.Select(
        geom_inter_pt, intermediate("sel_geom_pt", i, DEBUG), above
    )

    dissolve_fields = [in_fid_field, inter_cr_fid]

    valid_pt_path = intermediate("valid_stufe_pt", i, DEBUG)
    valid_pt = arcpy.management.Dissolve(
        sel_geom_pt, valid_pt_path, dissolve_fields, multi_part="MULTI_PART"
    )

    n_after_stufe = int(arcpy.management.GetCount(valid_pt).getOutput(0))
    print(f"  After STUFE filter: {n_after_stufe} geometry intersection points")

    if n_after_stufe == 0:
        print("  No candidates remain — skipping config")
        continue

    # Build: (input_oid, cross_oid) -> [point_geom, ...] for geometry intersections
    geom_inter_pts_lookup = {}  # (in_oid, cr_oid) -> [Point, ...]
    with arcpy.da.SearchCursor(valid_pt, [in_fid_field, inter_cr_fid, "SHAPE@"]) as cur:
        for in_oid, cr_oid, geom in cur:
            if geom is None:
                continue
            key = (in_oid, cr_oid)
            pts = []
            if geom.type == "multipoint":
                for p in range(geom.partCount):
                    pts.append(geom.getPart(p))
            else:
                pts.append(geom.firstPoint)
            geom_inter_pts_lookup.setdefault(key, []).extend(pts)

    # Valid (in_oid, cr_oid) pairs that passed STUFE rule
    valid_pairs = set(geom_inter_pts_lookup.keys())
    print(f"  Valid (input, crossing) pairs: {len(valid_pairs)}")

    # ---------------------------------------------------
    # 5. GENERATE SOLID MASKS — only for surviving features
    # ---------------------------------------------------

    print("  GENERATE SOLID MASKS (surviving features only)")

    # Select only the input features that appear in valid pairs
    valid_in_oids = set(k[0] for k in valid_pairs)
    valid_cr_oids = set(k[1] for k in valid_pairs)

    in_oid_str = ",".join(map(str, valid_in_oids))
    cr_oid_str = ",".join(map(str, valid_cr_oids))

    # Input masks
    input_mask_fc = f"{gdb}/QA_InputMasks_{i}"
    arcpy.management.SelectLayerByAttribute(
        input_aprx_lyr, "NEW_SELECTION", f"OBJECTID IN ({in_oid_str})"
    )
    original_input_cim = apply_solid_cim(input_aprx_lyr)
    aprx.save()
    try:
        arcpy.cartography.FeatureOutlineMasks(
            input_aprx_lyr,
            input_mask_fc,
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
    arcpy.management.SelectLayerByAttribute(input_aprx_lyr, "CLEAR_SELECTION")
    print(
        f"  Input masks: {int(arcpy.management.GetCount(input_mask_fc).getOutput(0))}"
    )

    # Crossing masks
    cross_mask_fc = f"{gdb}/QA_CrossMasks_{i}"
    arcpy.management.SelectLayerByAttribute(
        cross_aprx_lyr, "NEW_SELECTION", f"OBJECTID IN ({cr_oid_str})"
    )
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
    arcpy.management.SelectLayerByAttribute(cross_aprx_lyr, "CLEAR_SELECTION")
    print(
        f"  Cross masks: {int(arcpy.management.GetCount(cross_mask_fc).getOutput(0))}"
    )

    # ---------------------------------------------------
    # 6. INTERSECT MASKS (polygon + point)
    # ---------------------------------------------------

    print("  INTERSECT MASKS")

    # Confirm the FID field names in the mask FCs
    in_mask_fid = next(
        (
            f.name
            for f in arcpy.ListFields(input_mask_fc)
            if f.name.startswith("FID_DKM25")
        ),
        None,
    )
    cr_mask_fid = next(
        (
            f.name
            for f in arcpy.ListFields(cross_mask_fc)
            if f.name.startswith("FID_DKM25")
        ),
        None,
    )
    if not in_mask_fid or not cr_mask_fid:
        print(f"  WARNING: mask FID fields not found, skipping config {i}")
        continue

    mask_inter_poly = arcpy.analysis.Intersect(
        [input_mask_fc, cross_mask_fc],
        intermediate("mask_inter_poly", i, DEBUG),
        "ALL",
        output_type="INPUT",
    )
    mask_inter_pt = arcpy.analysis.Intersect(
        [input_mask_fc, cross_mask_fc],
        intermediate("mask_inter_pt", i, DEBUG),
        "ALL",
        output_type="POINT",
    )
    print(
        f"  Mask intersection points: {int(arcpy.management.GetCount(mask_inter_pt).getOutput(0))}"
    )

    # Confirm cr_fid in mask intersect output
    inter_cr_fid = next(
        (
            f.name
            for f in arcpy.ListFields(mask_inter_pt)
            if f.name.startswith("FID_DKM25") and f.name != in_mask_fid
        ),
        None,
    )
    if not inter_cr_fid:
        print(f"  WARNING: inter_cr_fid not found, skipping config {i}")
        continue

    # ---------------------------------------------------
    # 7. APPLY STUFE RULE ON MASK INTERSECT + DISSOLVE
    # ---------------------------------------------------

    print("  APPLY STUFE RULE ON MASK INTERSECT")
    sel_mask_pt = arcpy.analysis.Select(
        mask_inter_pt, intermediate("sel_mask_pt", i, DEBUG), above
    )
    sel_mask_poly = arcpy.analysis.Select(
        mask_inter_poly, intermediate("sel_mask_poly", i, DEBUG), above
    )
    print(
        f"  After STUFE filter: {int(arcpy.management.GetCount(sel_mask_pt).getOutput(0))}"
    )

    dissolve_fields = [in_mask_fid, inter_cr_fid]
    valid_pt = arcpy.management.Dissolve(
        sel_mask_pt,
        intermediate("valid_pt", i, DEBUG),
        dissolve_fields,
        multi_part="MULTI_PART",
    )
    valid_poly = arcpy.management.Dissolve(
        sel_mask_poly,
        intermediate("valid_poly", i, DEBUG),
        dissolve_fields,
        multi_part="MULTI_PART",
    )

    arcpy.management.AddField(valid_pt, "ORIG_OID", "LONG")
    with arcpy.da.UpdateCursor(valid_pt, ["OBJECTID", "ORIG_OID"]) as cur:
        for row in cur:
            row[1] = row[0]
            cur.updateRow(row)

    # ---------------------------------------------------
    # 8. BUILD LOOKUPS
    # ---------------------------------------------------
    mask_fid_to_ds_oid = {}
    with arcpy.da.SearchCursor(input_mask_fc, ["OBJECTID", in_mask_fid]) as cur:
        for mask_oid, src_fid in cur:
            mask_fid_to_ds_oid[mask_oid] = src_fid

    # Load mask geometries
    input_mask_geom = {}
    with arcpy.da.SearchCursor(input_mask_fc, [in_mask_fid, "SHAPE@"]) as cur:
        for fid, geom in cur:
            input_mask_geom[fid] = geom

    # Load line geometries
    # Replace the filtered geometry lookups with unfiltered ones

    input_line_geom = {}
    with arcpy.da.SearchCursor(input_ds, ["OBJECTID", "SHAPE@"]) as cur:
        for oid, geom in cur:
            input_line_geom[oid] = geom

    cross_line_geom = {}
    with arcpy.da.SearchCursor(cross_ds, ["OBJECTID", "SHAPE@"]) as cur:
        for oid, geom in cur:
            cross_line_geom[oid] = geom

    # Overlap polygon lookup: (in_mask_fid, cr_mask_fid) -> unioned polygon
    overlap_poly_lookup = {}
    with arcpy.da.SearchCursor(
        valid_poly, [in_mask_fid, inter_cr_fid, "SHAPE@"]
    ) as cur:
        for in_fid, cr_fid, poly_geom in cur:
            if cr_fid is not None and cr_fid != -1:
                key = (in_fid, cr_fid)
                overlap_poly_lookup[key] = (
                    overlap_poly_lookup[key].union(poly_geom)
                    if key in overlap_poly_lookup
                    else poly_geom
                )

    # Point count lookup from dissolved multipoint
    pt_count_lookup = {}
    with arcpy.da.SearchCursor(valid_pt, [in_mask_fid, inter_cr_fid, "SHAPE@"]) as cur:
        for in_fid, cr_fid, geom in cur:
            if cr_fid is not None and cr_fid != -1 and geom is not None:
                pt_count_lookup[(in_fid, cr_fid)] = geom.partCount

    # ---------------------------------------------------
    # 9. DETECT TRUE CROSSINGS
    # ---------------------------------------------------

    print("  DETECT TRUE CROSSINGS")

    endpoint_oids = set()
    no_polygon_oids = set()
    true_cross_oids = set()

    count_no_poly = 0
    count_too_few_pts = 0
    count_no_geom_inter = 0
    count_true = 0
    count_endpoint = 0

    with arcpy.da.SearchCursor(
        valid_pt, ["ORIG_OID", in_mask_fid, inter_cr_fid]
    ) as cur:
        for orig_oid, in_fid, cross_fid in cur:
            if cross_fid is None or cross_fid == -1:
                continue

            key = (in_fid, cross_fid)
            overlap_poly = overlap_poly_lookup.get(key)

            # RULE 1 — no polygon
            if not overlap_poly:
                count_no_poly += 1
                no_polygon_oids.add(orig_oid)
                continue

            # RULE 2 — fewer than 4 mask boundary points
            if pt_count_lookup.get(key, 0) < 4:
                count_too_few_pts += 1
                endpoint_oids.add(orig_oid)
                continue

            # RULE 3 — geometry intersect pre-filter:
            # at least one geometry crossing point must exist for this pair
            geom_pts = geom_inter_pts_lookup.get((in_fid, cross_fid), [])
            if not geom_pts:
                count_no_geom_inter += 1
                endpoint_oids.add(orig_oid)
                continue

            # RULE 4 — per-part erase + endpoint + geometry containment check
            in_mask = input_mask_geom.get(in_fid)
            in_line = input_line_geom.get(in_fid)
            cr_line = cross_line_geom.get(cross_fid)

            crossing, valid_parts = (
                is_true_crossing(in_mask, in_line, cr_line, overlap_poly, geom_pts)
                if in_mask and in_line and cr_line
                else (False, [])
            )

            if crossing and valid_parts:
                overlap_poly_lookup[key] = keep_parts(overlap_poly, valid_parts)
                count_true += 1
                true_cross_oids.add(orig_oid)
            else:
                count_endpoint += 1
                endpoint_oids.add(orig_oid)

    print(f"  No polygon:              {count_no_poly}")
    print(f"  Too few pts (<4):        {count_too_few_pts}")
    print(f"  No geometry intersect:   {count_no_geom_inter}")
    print(f"  True crossings:          {count_true}")
    print(f"  Endpoint touches:        {count_endpoint}")

    # ---------------------------------------------------
    # 10. BUILD OUTPUT
    # ---------------------------------------------------

    print(f"  BUILD OUTPUT")
    out_fc_i = f"{gdb}/QA_ValidCrossings_{i}"

    if not arcpy.Exists(out_fc_i):
        arcpy.management.CopyFeatures(valid_pt, out_fc_i)
        for fname, ftype, flen in [
            ("CROSS_SOURCE", "TEXT", 30),
            ("CROSS_TYPE", "TEXT", 20),
            ("MASK_ERROR", "SHORT", 0),
            ("MASK_ERROR_TYP", "TEXT", 50),
            ("PART_COUNT", "SHORT", 0),
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
        in_mask_fid,
        inter_cr_fid,
        "CROSS_SOURCE",
        "CROSS_TYPE",
        "MASK_ERROR",
        "MASK_ERROR_TYP",
        "PART_COUNT",
    ]

    with arcpy.da.UpdateCursor(out_fc_i, update_fields) as cur:
        for row in cur:
            shape = row[0]
            oid = row[1]
            line_fid = row[2]
            cross_fid = row[3]
            row[4] = cross_ds
            row[8] = shape.partCount if shape else 0

            if oid in no_polygon_oids:
                row[5] = "NO_POLYGON"
                row[6] = 0
                row[7] = "N/A"
            elif oid in endpoint_oids:
                row[5] = "ENDPOINT_TOUCH"
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
