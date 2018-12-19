#!/bin/bash
#aa=data
#aa=A
aa=1
b=0.25
#c=ChallengeA
c=LightConeHector
#c=ChallengeJapanCMASS2

d=0
#e=ChallengeHDvFFT_IR06
#e=PTchallengeCMASS2
e=LightConeHectorPatchyWideHDvFFT_IR06
#e=LightConeHectorDataHDvFFT_IR06
#f=0
f=0.08
g=1
h=0
ii=100

cd /exports/pierre/EFTofBOSS

python2 chi2.py $aa $b $c $d $e $f $g $h $ii
#python2 verifchi2.py $aa $b $c $d $e $f $g $h $ii
#python2 verifchi2fidcosmo.py $aa $b $c $d $e $f $g $h $ii

wait

