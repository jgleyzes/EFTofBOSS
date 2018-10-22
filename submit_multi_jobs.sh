#submit z=0.61 jobs first, simtype is LightConeDida
#Simnames 

#simname=ChallengeA box name is conChallengeAwideh
#simname=LightConeHector box name is conLightConeHectorwideh
planck=0
simname=LightConeHector
marg=1
theoryCovness=0
#for kmax in 0.1
for kmax in 0
#for kmax in 0.25
do
#for kmaxbisp in 0.07 0.1
for kmaxbisp in 0.07
do
#for number in A B F G 
for number in 5 8  
do
 #the above is for challenge boxes
sbatch --job-name=b$number.$kmax.$simname --output=outfiles/b$number.$kmax.$kmaxbisp.$simname.Planck$planck.out cluster_MCMC.sh $number $kmax $simname $planck conLightConeHectorwideh $kmaxbisp $marg $theoryCovness
done
#sbatch --job-name=bdata.$kmax.$simname --output=outfiles/bdata.$kmax.$kmaxbisp.$simname.Planck$planck.out cluster_MCMC.sh data $kmax $simname $planck conLightConeHectorwideh $kmaxbisp $marg
done
done

