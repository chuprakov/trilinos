#!/bin/sh -f

if [ -a $1 ]; then

		echo 'no need to make directory ' $1

else

		echo 'making directory ' $1
		/bin/mkdir -p $1

fi




