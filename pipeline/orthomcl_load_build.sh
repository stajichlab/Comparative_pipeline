

module load orthomcl

orthomclLoadBlast dbconfig.txt Dothideomycetes.BLASTP.bpo
orthomclPairs dbconfig.txt  pairs.log cleanup=yes
orthomclDumpPairsFiles dbconfig.txt

mcl mclInput --abc -I 1.5 -o Dothideomycetes.OrthoMCL.out

# orthomclMclToGroups
