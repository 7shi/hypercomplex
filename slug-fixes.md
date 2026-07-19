# 旧slug修正チェックリスト

`uv run slug_fixes.py` で生成。リネーム表は `slug_fixes.py`、決定の経緯は `SLUG.md` を参照。

記事ごとに、Mathlog側の**本文**（`[[slug]]` の置換）と**参考文献登録**（ラベルの付け直し）、およびローカルの md と `refs/{ID}.toml` を修正する。

## clif/01-representation.md — Mathlog未公開（ローカルのみ）

- [ ] `pauli-qua` → `7shi-bq`（本文3箇所）
- [ ] `qua-tensor` → `7shi-qt`（本文5箇所）

## clif/03-gpm-gca.md — Mathlog未公開（ローカルのみ）

- [ ] `clif-rep` → `7shi-clif1`（本文7箇所）
- [ ] `clif-sc` → `7shi-clif2`（本文1箇所）
- [ ] `gca` → `wiki-gca`（本文2箇所）
- [ ] `gpm` → `wiki-gpm`（本文2箇所）
- [ ] `pauli-qua` → `7shi-bq`（本文3箇所）
- [ ] `qua-nonion` → `7shi-nonion`（本文2箇所）
- [ ] `qua-tensor` → `7shi-qt`（本文3箇所）

## clif/04-weyl-algebra.md — Mathlog未公開（ローカルのみ）

- [ ] `gca` → `wiki-gca`（本文2箇所）
- [ ] `gpm` → `wiki-gpm`（本文1箇所）
- [ ] `lie-u1` → `7shi-lie1`（本文1箇所）
- [ ] `qua-nonion` → `7shi-nonion`（本文1箇所）

## clif/05-schrodinger.md — Mathlog未公開（ローカルのみ）

- [ ] `gca` → `wiki-gca`（本文2箇所）
- [ ] `gpm` → `wiki-gpm`（本文3箇所）
- [ ] `nelson-sqm` → `7shi-leb6`（本文1箇所）
- [ ] `weyl-alg` → `7shi-clif4`（本文2箇所）

## hopf/01-quaternion.md — https://mathlog.info/articles/enXfjr4Xvas6pfXLQddc

- [ ] `7shi` → `7shi-qrot`（本文1箇所、参考文献登録）
- [ ] `7shi-hc` → `7shi-homog`（本文1箇所、参考文献登録）
- [ ] `wikipedia` → `wiki-hopf`（本文4箇所、参考文献登録）

## hopf/02-spinor-tensor.md — https://mathlog.info/articles/co6jC5vshPXR0EgmAjX3

- [ ] `7shi-c` → `7shi-coord`（本文1箇所、参考文献登録）
- [ ] `wikipedia-dp` → `wiki-dp`（本文1箇所、参考文献登録）

## hopf/04-extension.md — Mathlog未公開（ローカルのみ）

- [ ] `7shi-oct` → `7shi-oct1`（本文1箇所）

## hopf/05-entanglement.md — Mathlog未公開（ローカルのみ）

- [ ] `7shi-bd` → `7shi-bloch`（本文3箇所）
- [ ] `7shi-ext` → `7shi-hopfext`（本文6箇所）
- [ ] `7shi-oct` → `7shi-oct1`（本文1箇所）
- [ ] `7shi-tensor` → `7shi-tp`（本文1箇所）
- [ ] `lie-su3` → `7shi-lie5`（本文4箇所）

## hopf/06-clifford-gates.md — Mathlog未公開（ローカルのみ）

- [ ] `7shi-bd` → `7shi-bloch`（本文9箇所）
- [ ] `7shi-ent` → `7shi-entangle`（本文5箇所）
- [ ] `7shi-rot` → `7shi-clrt`（本文1箇所）
- [ ] `7shi-spin` → `7shi-lie3`（本文1箇所）
- [ ] `7shi-su2` → `7shi-lie2`（本文4箇所）
- [ ] `7shi-tensor` → `7shi-tp`（本文1箇所）
- [ ] `gpm` → `wiki-gpm`（本文1箇所）
- [ ] `pauli-qua` → `7shi-bq`（本文4箇所）
- [ ] `weyl-alg` → `7shi-clif4`（本文1箇所）

## hopf/c2-to-s2.md — https://mathlog.info/articles/OHZI3f5bQ4hB0S7dn6lv

- [ ] `7shi-hopf1` → `7shi-h`（本文1箇所、参考文献登録）

## hopf/homogeneous.md — https://mathlog.info/articles/f5v2F6yOFhqtnpDtFtcP

- [ ] `7shi-qhopf1` → `7shi-h`（本文1箇所、参考文献登録）

## lie/02-su2-so3.md — https://mathlog.info/articles/Utdur1fLLzrWVHOJHifj

- [ ] `7shi-qmat` → `7shi-qcm`（本文1箇所、参考文献登録）

## lie/03-spin.md — Mathlog未公開（ローカルのみ）

- [ ] `7shi-hopf1` → `7shi-h`（本文1箇所）
- [ ] `7shi-su2` → `7shi-lie2`（本文1箇所）

## lie/04-spin4.md — Mathlog未公開（ローカルのみ）

- [ ] `7shi-spin` → `7shi-lie3`（本文3箇所）
- [ ] `7shi-su2` → `7shi-lie2`（本文1箇所）

## misc/spherical-coords.md — https://mathlog.info/articles/2HkSkfawcYJYgE6rVgu6

- [ ] `7shi` → `7shi-h`（本文1箇所、参考文献登録）
- [ ] `7shi-2` → `7shi-s`（本文1箇所、参考文献登録）
- [ ] `wikipedia` → `wiki-nsphere`（本文1箇所、参考文献登録）
- [ ] `wikipedia-s3` → `wiki-s3`（本文1箇所、参考文献登録）

## oct/01-octonion.md — https://mathlog.info/articles/vizAwdiKusmQ6pdsPTLE

- [ ] `oct-7rot` → `7shi-7rot`（本文1箇所、参考文献登録）
- [ ] `oct-cl6` → `7shi-cl6`（本文1箇所、参考文献登録）
- [ ] `oct-nonassoc` → `7shi-nonassoc`（本文1箇所、参考文献登録）

## oct/02-7d-3rot.md — https://mathlog.info/articles/cWj1bUNwa7E3kgMhvJb0

- [ ] `7shi-rot` → `7shi-clrt`（本文1箇所、参考文献登録）
- [ ] `oct` → `7shi-oct1`（本文1箇所、参考文献登録）

## oct/03-oct-left-mul.md — https://mathlog.info/articles/PZJuaRWHYidsIg9u83j8

- [ ] `oct` → `7shi-oct1`（本文1箇所、参考文献登録）
- [ ] `ref1` → `hasebe`（本文1箇所、参考文献登録）
- [ ] `ref2` → `tian`（本文1箇所、参考文献登録）

## oct/nonassociativity.md — https://mathlog.info/articles/fXYD8gdyHdnyPvoFnWPN

- [ ] `7shi-oct480` → `7shi-480`（本文3箇所、参考文献登録）

## qua/01-pauli-qua.md — https://mathlog.info/articles/lZ1X3t6exNS3NrNArqji

- [ ] `7shi-hopf` → `7shi-h`（参考文献登録）
- [ ] `7shi-pcoord` → `7shi-coord`（本文1箇所、参考文献登録）

## qua/02-nonion.md — https://mathlog.info/articles/3gMUAXJD714J2LtlfIqU

- [ ] `gca` → `wiki-gca`（本文1箇所、参考文献登録）
- [ ] `gpm` → `wiki-gpm`（本文1箇所、参考文献登録）

## qua/cd/grok45.md — Mathlog未公開（ローカルのみ）

- [ ] `7shi-su2` → `7shi-lie2`（本文3箇所）

## qua/real-pauli-trace.md — https://mathlog.info/articles/Ale2XUWFasqCiU2V5Q7p

- [ ] `7shi-pbq` → `7shi-bq`（本文1箇所、参考文献登録）

## qua/spinor-ideal.md — Mathlog未公開（ローカルのみ）

- [ ] `7shi-oct-left-mul` → `7shi-cl6`（本文3箇所）
- [ ] `clif-rep` → `7shi-clif1`（本文2箇所）

## qua/tensor-from-complex.md — https://mathlog.info/articles/OunyiPffmIKhwwdmSrsA

- [ ] `7shi-pdq` → `7shi-bq`（本文1箇所、参考文献登録）

## vec-oct/01-reflection.md — https://mathlog.info/articles/yZGOwcB2XXpfEvlCuSW4

- [ ] `7shi-ort` → `7shi-7rot`（本文1箇所、参考文献登録）

## vec-oct/geometric-product-exp.md — https://mathlog.info/articles/YZmxak6ObeP6rLQnyU2V

- [ ] `7shi-pdq` → `7shi-bq`（本文1箇所、参考文献登録）
- [ ] `7shi-vo` → `7shi-mir`（本文1箇所、参考文献登録）

