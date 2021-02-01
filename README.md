# mcr_with_bemg
pyMCR with BEMG constraint implimentaion example
pyMCRをつかってLC-PDA (LC-DAD)クロマトピークを分離する際に、BEMGモデル関数制約を入れるサンプルコード

 - do_mcr.py : the MCR example　　MCRを実行します
 - constraint.py:  BEMG constraint 　BEMG制約
 - make_dummy_lc_pda.py: make dummy data for experiment　　実験用のダミーデータを生成します
 - sim_chrom.py : chromatogram simulation code 　クロマト波形をシミュレーションします
 - isotherm.py : isotherm for chromatogram simulation　クロマト波形シミュレーション用の吸着等温線



pyMCRは初期値が必要なのでFastICAで初期値を作る
![Fig1](https://github.com/akirayou/mcr_with_bemg/blob/main/img/Figure_1.png)

そのままpyMCRを実行するとそれなりには分離してくれる
![Fig2](https://github.com/akirayou/mcr_with_bemg/blob/main/img/Figure_2.png)

BEMG制約を入れると、クロマト形状がBEMG関数をとる形のモノに絞ってくれる。

![Fig3](https://github.com/akirayou/mcr_with_bemg/blob/main/img/Figure_3.png)

クロマトピーク形状をシミュレーションで作っているので、そこそこリアルなリーディング・テーリングを楽しめます。
実装サンプルでは段理論(Plate therory)によるクロマトシミュレーションをRadke-Prausnitz等温線を使って、
実行していますが、FreundlichやLangmuirも選べます。また配管デッドボリュームによる一次遅れテーリングも付加できます。

#Refferences 参考文献
 - https://pages.nist.gov/pyMCR/
 - https://www.sciencedirect.com/science/article/abs/pii/S0021967316312535
 - https://patents.google.com/patent/WO2016035167A1/ja
