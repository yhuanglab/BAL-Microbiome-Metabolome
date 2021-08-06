#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=DIABLO
#SBATCH --mail-user=madapoos@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000m 
#SBATCH --time=48:00:00
#SBATCH --account=yvjhuang99
#SBATCH --partition=standard

# The application(s) to execute along with its input arguments and options:

module load R

for i in /home/madapoos/Spiromics_No_Severe_No_Healthy_Reg/Output/PEX_HCUDRUG_CAT_365.y_1_Oneormore_0_None	/home/madapoos/Spiromics_No_Severe_No_Healthy_Reg/Output/FVC_BDRESPONSE_PCT_ctb_1_12above
do
		echo "***********"
		echo "Starting $i"
		cd $i #Enter folder
		Rscript /home/madapoos/Spiromics_No_Severe_No_Healthy_Reg/Code/3_DIABLO_Cluster_Version.R
		echo "Modeling Completed"
		Rscript /home/madapoos/Spiromics_No_Severe_No_Healthy_Reg/Code/3_DIABLO_Feature_Aggregation.R
		echo "$i Completed"
		echo "***********"
		cd .. #Exit folder
done

/bin/hostname
sleep 60