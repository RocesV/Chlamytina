
#!bin/bash

for STATE in {2..40..1}
do
   mkdir ../ChromHMM_Outputs/$STATE
   for ITER in {1..11..1}
   do
     mkdir ../ChromHMM_Outputs/$STATE/$ITER
     echo "Running $STATE model, $ITER iteration | $(date)"
     if [ $ITER -eq 1 ]; then 
       java -mx20000M -jar ChromHMM.jar LearnModel -b 150 -init information ../binariesV/ ../ChromHMM_Outputs/$STATE/$ITER/ $STATE Cre55
     else
       java -mx20000M -jar ChromHMM.jar LearnModel -b 150 -init random ../binariesV/ ../ChromHMM_Outputs/$STATE/$ITER/ $STATE Cre55
     fi
   done		 	
done 


