# OUMasken_QA 
## Overview

**Cartographic Crossing Mask Detection**

This tool detects locations where a road or railway passes beneath another road or railway at a different vertical level (`STUFE`) and checks whether an appropriate cartographic mask polygon exists at each crossing.

The purpose is to identify missing or misplaced masks that can lead to cartographic conflicts in the map.

The repository contains the original crossing QA implementation and a geometry-based test implementation developed to investigate an alternative processing approach.

## Repository Contents

| File | Description |
|---|---|
| `Crossing_QA.py` | Original implementation. Generates and intersects feature masks to identify potential crossings and check existing cartographic masks. |
| `Crossing_QA_geom_test.py` | Experimental implementation that identifies candidate crossings using centre-line geometry intersections before generating masks. |
| `Crossing_detection.pptx` | Presentation slides explaining the crossing detection approach. |
| `DKM25_CartoProcess_Config.cp.xml` | Cartographic process configuration used to test the scripts. |

### Test data

The scripts were tested using the **`WU1207_Thun` work unit**.

The associated geodatabase (GDB) and ArcGIS Pro project (`.aprx`) are not included in this repository because of their size. 
The scripts therefore require access to the corresponding work-unit data and a compatible ArcGIS environment.

## Processing Configuration

Both implementations use the cartographic process configuration defined in `DKM25_CartoProcess_Config.cp.xml`.

Each relevant process configuration defines the input features, crossing features, vertical-level rule, and existing cartographic masks used for QA.

### Inputs per process configuration

- **Input dataset and WHERE clause:** the road or railway being evaluated as the lower feature.
- **Crossing dataset(s) and WHERE clause:** the road or railway potentially passing above the input feature.
- **`AboveExpression`:** SQL WHERE clause defining the vertical-level relationship, typically `STUFE_input < STUFE_crossing`.
- **`MaskDataset`:** feature class containing existing cartographic mask polygons.
- **`MaskAssociation`:** table linking road features (`LINEID`) to their masks (`MASKID`).

The configuration determines which feature combinations are evaluated and which vertical-level relationships qualify as potential crossings.

## Pipeline: `Crossing_QA.py`

The original implementation uses mask generation and polygon intersections to detect crossings.

1. **Select inputs and generate input masks**
2. **Load QA mask tables**
3. **Process each configuration**
   1. Select crossings and generate crossing masks.
   2. Intersect the masks.
   3. Apply the `AboveExpression` rule and dissolve the results.
   4. Build lookup tables.
   5. Filter for true crossings.
   6. Build the output.
4. **Merge results and generate a summary**

This approach uses generated mask geometry as an important part of the crossing detection process.

## Pipeline: `Crossing_QA_geom_test.py`

The geometry-based implementation was developed as an alternative to the original approach. It uses centre-line intersections to identify candidate crossing pairs before generating masks, with the aim of reducing unnecessary mask generation.

1. **Select input roads and crossing roads/railways** from the geodatabase.
2. **Find candidate crossing pairs** using centre-line geometry intersections.
3. **Apply the `STUFE` rule**, retaining only pairs where the input feature is at a strictly lower level than the crossing feature.
4. **Generate solid Feature Outline Masks** for the surviving features. Dash effects are temporarily removed from the Cartographic Information Model (CIM) during mask generation.
5. **Intersect the generated masks** to identify visual overlap areas.
6. **Filter for true crossings** using the three conditions implemented in `is_true_crossing`.
7. **Check existing cartographic masks** and flag missing or misplaced masks.

By filtering candidate pairs before mask generation, this approach aims to avoid generating masks for features that cannot form a valid crossing under the configured rules.

## Outputs

The scripts generate QA point feature classes in the geodatabase.

| Output | Description |
|---|---|
| `QA_ValidCrossings_<i>` | Per-configuration output containing one row per detected crossing pair. |
| `QA_ValidCrossings` | Merged result across all processed configurations. |

### Output fields

| Field | Description |
|---|---|
| `CROSS_SOURCE` | Identifies the crossing dataset that produced the result. |
| `CROSS_TYPE` | Classification of the detected crossing: `CROSSING`, `ENDPOINT_TOUCH`, or `NO_POLYGON`. |
| `MASK_ERROR` | Mask QA status: `0` = OK, `1` = error. |
| `MASK_ERROR_TYP` | Mask QA classification: `OK`, `NO_MASK`, `MASK_MISPLACED`, or `N/A`. |

The output can be used to locate potential cartographic mask issues and compare detection results between implementations.

## Comparison of the Two Implementations

The two scripts do not currently produce identical results.

In general, the geometry-based implementation (`Crossing_QA_geom_test.py`) misses some crossing cases identified by the original implementation. However, in configuration 2 **`Create CrossingMasks STRASSE einlinig unter STRASSE, EISENBAHN einlinig (Dec)`** the geometry-based implementation identifies substantially more crossings.

These differences have not yet been fully investigated yet.
