
# Need to be on uvcontsub_exclude_freqrange for now
git clone https://github.com/e-koch/phangs_imaging_scripts.git
git checkout origin/uvcontsub_exclude_freqrange
git switch -c uvcontsub_exclude_freqrange


/opt/casa-6.2.1-7-pipeline-2021.2.0.128/bin/casa --nogui --log2term -c ~/lglbs_hi_scripts/run_lglbs_staging_dwarfctr.py ic10ctr
/opt/casa-6.2.1-7-pipeline-2021.2.0.128/bin/casa --nogui --log2term -c ~/lglbs_hi_scripts/run_lglbs_staging_dwarfctr.py ic1613ctr
/opt/casa-6.2.1-7-pipeline-2021.2.0.128/bin/casa --nogui --log2term -c ~/lglbs_hi_scripts/run_lglbs_staging_dwarfctr.py wlmctr
# ngc6822 has no centre only pointing


