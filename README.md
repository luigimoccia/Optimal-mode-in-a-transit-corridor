# Optimal-mode-in-a-transit-corridor

The python codes in this directory solve models for the optimization of a transit line 
and for the strategic assessment of mode choice presented in: 

L. Moccia and G. Laporte, “Improved models for technology choice in a transit corridor with fixed demand”, 2015

A draft of this paper is in this directory as “paper.pdf”


COMMENTS
Please write to: 
moccia@icar.cnr.it


REQUIREMENTS
 Python 2.7 with Numpy, Scipy, Shapely, Matplotlib

HOW TO USE
To solve model I type from the command line:
python code_mod1.py

Similarly for the other models, e.g. 
python code_mod2.py

Etc.

INPUT
The input parameters are in the file common_parameters.py.
See the comments on this file for details.

OUTPUT
The program reports some self-explanatory computational details at screen,
and yields several figures as pdf files.
For example, code_mod1.py generates pdf files named 'fig_M1_*.pdf'
On a mac platform these files are opened upon completion of the code.



