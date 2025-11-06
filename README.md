
# Read truncated .acq files
The script opens .acq files using the Bioread package (https://pypi.org/project/bioread/). When power is cut or due to technical reasons the data collection stops mid-read, event markers are not written to the file (this takes place after data finalization by AckNowledge). However, since the data is streamed from the BIOPAC data collection system, the data is still saved. Through reading the information from the channel headers we get from Bioread, we can manually read and cast the data. Finally, the data is exported as a .csv file for further processing/analysis.

![alt text](https://github.com/freesynapse/corrupt_bioread/blob/main/img/sanity.png)