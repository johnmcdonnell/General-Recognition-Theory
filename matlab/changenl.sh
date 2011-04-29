#!/bin/sh

# Just runs the endline fixing vim script on all the matlab files.

for file in *.m; 
do 
	vi -s changenl.vim $file; 
done

for file in GRTdemos/*.m; 
do 
	vi -s changenl.vim $file; 
done
