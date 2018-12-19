#Simnames 
#simname=LightConeHector box name is conLightConeHectorFFT
planck=0
marg=1
theoryCovness=0.008
simname=LightConeHector
grid=LightConeHectorFFTwideh
#simname=ChallengeA
#grid=conChallengeAwideh

numnode=0
prior=100

cd /exports/pierre/EFTofBOSS


for kmax in 0.25 0.3 0.35 0.4 0.45 0.5
do
#for kmaxbisp in 0.07 0.1
for kmaxbisp in 0
do
for number in 8
#for number in {1..9} 
do
numnode=$((numnode+1))
 #the above is for challenge boxes
qsub -l nodes=node08:ppn=8 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
#qsub -k oe -l nodes=node10:ppn=48 -N bdata.$kmax.$simname -o log/output_bdata.$kmax.$kmaxbisp.$simname.Planck$planck.log -e log/error_bdata.$kmax.$kmaxbisp.$simname.Planck$planck.log -v aa=data,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg cluster_MCMC_Mercury.sh
done

done
