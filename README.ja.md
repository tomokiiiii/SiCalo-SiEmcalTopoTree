# SiCalo documentation

この directory の simulation と `SiEmcalTopoTree` 作成手順は、入力 sample ごとに分けています。

- [Pythia pp simulation からの作成](README_PYTHIA.ja.md)
- [single-particle gun simulation からの作成](README_SINGLEGUN.ja.md)

このrepositoryでは、MLに渡すROOT fileが作られるまでを説明します。
ElectronClassifier本体と学習・評価・filter・plotのコードは公開対象に含めていません。

## 現在の使い分け

| 目的 | 使用する入口 |
|---|---|
| 現行 Pythia position v5 production | `dstmacro/sim/run.sh` |
| ElectronClassifier で使用中の pp v4 sample を再現 | `Fun4All_PYTHIA_SiliconRecov2.C` |
| electron single gun | `macros/run_singlev2.sh` |
| negative-pion single gun | `macros/run_single.sh` |
| particle を指定した single gun | `Fun4All_singleParticle_SiliconPositionv4.C` を直接実行 |

新しい macro では、一回の実行で Fun4All analysis、`ana_*.root`、`make_SiEmcalTopoTree`、ML feature の追加まで連続して行います。世代による例外と output 名は各文書に記載しています。

