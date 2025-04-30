#!/bin/bash
source activate base

# Run the python script with different arguments in parallel using nohup
nohup python ./AFD_Word_Scrape.py WR > output_WR.log 2>&1 &
nohup python ./AFD_Word_Scrape.py CR > output_CR.log 2>&1 &
nohup python ./AFD_Word_Scrape.py ER > output_ER.log 2>&1 &
nohup python ./AFD_Word_Scrape.py SR > output_SR.log 2>&1 &

# Wait for all background processes to finish
wait

echo "All scripts have finished running."