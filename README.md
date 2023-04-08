# gradient_method_triangles
最急降下法を用いて2つの三角形の最大重ね合わせを求める  
（プログラミング課題２)


---
## Contents
File/Folder | Description
-|-
``task2_0221.c`` | ``input1.txt``に記載された2つの三角形の座標を読み込み、重ね合わせが最大になるような配置の座標および画像ファイルをoutputする
``input1.txt`` | ２つの三角形の座標が書かれたinputファイル


---
## Description
- ### task2_0221
 ２つの三角形の各辺の中点を合わせた計９パターンを初期値とし、x軸, y軸, 角度の３変数を用いて最急降下法で局所的最適解を導出。局所的最適解同士を比較し、最も重ね合わさった面積が大きい解を最適解として``output3.txt``および``output.png``で出力する。  


Outline  
1. inpitファイルから座標を読み込む
1. 初期値に移動
1. 数値微分から傾きを出し、重なる面積が大きくなる方向に移動させる
1. 重なった面積の最大値が更新されなくなるまで続ける
1. 最大値の更新が止まった値を局所的最適解として再び新たな初期値から最急降下法を行う
1. それぞれの局所的最適解を比較して最も重なった面積が大きいものを最適解として出力する  


（なお、重なった面積は、頂点と交点の集合から相手三角形の内側に存在するものを抽出し偏角ソートを行ったうえで座標法を用いて導出する。）


- ### input.txt
```
(0.0,0.0)(7.0,4.0)(9.0,0.0)
(0.0,0.0)(2.0,4.0)(6.0,0.0)
```
２つの三角形の座標を(x座標,y座標)のように記入する  
記入する座標は時計周りの順番で記入する。（順番が間違っていると実行時に``INPUT ERROR``と表示される）
