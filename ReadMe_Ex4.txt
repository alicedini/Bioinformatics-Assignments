Dini Alice - Assignment no. 4

The two files Ex4.py and Ex4_models.py provide possible implementations to:
- Compute the inside and outside CpG island models 
	based on CpG islands in chromosome 22 (which are given in attachment as well,
	previously computated with Ex4_models.py)
- Score a query sequence by means of those models
- Scan a long genome sequence using a sliding window, to obtain which 
	of these windows in the query genome have a high probability of 
	being inside a CpG island.

The file which has to be run is Ex4.py (possibly with Python 3.7.x): it has 
to be stored together with Ex4_models.py and the two files containing the models, 
"inside_file.npy" and "outside_file.npy".
To run it, it is enough to locate yourself into the directory where all the files 
are stored using the shell, and execute "python Ex4.py". 
The program is completely case insensitive. 
If the user wants to recompute the models, it is enough to delete the .npy file provided,
or to run models() function from Ex4_models.py; in this latter case, the .npy models will
be substituted and recomputed.

The user has the possibility to choose through the user interface either to score a 
query sequence, or to scan its query using a sliding window : in this latter case, 
if the input sequence is smaller than the average CpG island length, the user itself has 
the possibility to decide the size of the sliding window. Since the closer it is the length
of the query genome to the average CpG length, the smaller the number of CpG islands considered,
to obtain interesting results a considerable number of nucleotides should make up the input.

Case 1: Scoring a query sequence
The output will be the score of the sequence given the inside and the outside models, 
along with a prediction about being either inside or outside a CpG island.

Case 2: Scanning a long genome sequence using a sliding window
The output will be a tabulated list of offsets and scores referring only to those
regions with size equal to the window that returned a positive score, meaning they 
would be probably found inside a CpG island. 

Ex4_models.py provides the implementation of the code required to build the inside and 
the outside models. They are going to be saved in .npy format, and automatically read by 
Ex4.py when required; if they are not found in the directory where Ex4.py is stored, 
the models are recomputed.
Several libraries are required for the proper functioning of the module, namely:
- numpy
- random
- statistics
- re
- gzip
- urllib.

numpy library has been exploited for its methods devoted
to arrays' uploading and downloading from the file itself, and generally to deal with 
matrices in an optimal way. 

random library is used to generate random points in chromosome 22 to obtain DNA sequences for 
the computation of the outside model.

statistics library has been used during the coding session to compute the average CpG length.

re library, for regular expression, is used to compute the frequencies of each dinucleotide,
since re.findall method, differently from .count takes into account overlaps.

urllib is used to fetch the files using remote access.
Since the program directly fetches the files it requires directly from the servers where
they are stored, an internet connection is also required.
In fact, this first modules proceeds by opening the url where the files are available,
unzipping them in the case of chromosome 22 FASTA file, and decoding them into UTF-8 format,
since they are downloaded as strings of bytes.
A more detailed description about what the code actually does is available in the 
module in the commented parts. 

Ex4.py provides the user interface and the accomplishment of the two cases exposed above,
hence given a sequence s, it computes its score according to the models produced with
Ex4_models.py, by loading the .npy files created. The same occurs when of that sequence s,
the user wants a complete scan of s using a sliding window, whose size can vary if s length is
less than 567 (average CpG length).
Also here several libraries are used:
- numpy
- math
- os
- tabulate.

numpy is here used for the same reasons as Ex4_models.py's ones.

math is used to compute the logarithm of the ratio between the probability of being inside
and the probability of being outside a CpG island

os is used to check if the models exist in the directory where Ex4.py is found; if not so,
the modules get de novo computed.

tabulate is used to provide a user-friendly layout for the output of the offsets and scores
in Case 2.
A more detailed description about what the code actually does is available also in this
case in the module in the commented parts. 



