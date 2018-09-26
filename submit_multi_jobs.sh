#submit z=0.61 jobs first, simtype is LightConeDida
#Simnames 

#simname=ChallengeA
#should chain have planck rd prior?
planck=1
simname=LightConeHector

#for kmax in 0.1
for kmax in 0.25
#for kmax in 0.25
do
for kmaxbisp in 0.07 0.09 
#for kmaxbisp in 0.07
do
#for number in A B F 
for number in {1..9}
do
 #the above is for challenge boxes
sbatch --job-name=b$number.$kmax.$simname --output=outfiles/b$number.$kmax.$kmaxbisp.$simname.Planck$planck.out cluster_MCMC.sh $number $kmax $simname $planck conLightConeHectorwideh $kmaxbisp
done
sbatch --job-name=bdata.$kmax.$simname --output=outfiles/bdata.$kmax.$kmaxbisp.$simname.Planck$planck.out cluster_MCMC.sh data $kmax $simname $planck conLightConeHectorwideh $kmaxbisp
done
done

