#!/bin/sh
#
#  Run tests for TSF
#         
#################################################################################

#################################################################################
#
# Some functions
#
#################################################################################

write_html_header() {
echo "<!DOCTYPE HTML PUBLIC \"TSF Test Results\">" >> $TSF_TEST/tsftest.html
echo "<HTML> " >> $TSF_TEST/tsftest.html
echo " <HEAD> " >> $TSF_TEST/tsftest.html
echo "	 <TITLE>TSF Test Results</TITLE> " >> $TSF_TEST/tsftest.html
echo " </HEAD> " >> $TSF_TEST/tsftest.html
echo " <BODY> " >> $TSF_TEST/tsftest.html
return 0;
}

write_table_header() {
echo "" >> $TSF_TEST/tsftest.html
echo "   <P><B><FONT SIZE="+2">Results from TSF test</FONT></B> </P> " >> $TSF_TEST/tsftest.html
echo "   <P>Run on: $1</P> " >> $TSF_TEST/tsftest.html
shift
# echo "	 <P>Users who made commits: $*</P> " >> $TSF_TEST/tsftest.html
echo "   <TABLE BORDER="1"> " >> $TSF_TEST/tsftest.html
echo "     <TR> " >> $TSF_TEST/tsftest.html
echo "       <TD><FONT SIZE="+2"><B>Test Name</B></FONT></TD> " >> $TSF_TEST/tsftest.html
echo "       <TD><FONT SIZE="+2"><B>Status</B></FONT></TD> " >> $TSF_TEST/tsftest.html
echo "       <TD><FONT SIZE="+2"><B>Comments</B></FONT></TD> " >> $TSF_TEST/tsftest.html
echo "     </TR> " >> $TSF_TEST/tsftest.html
return 0;
}

write_table_footer() {
echo "   </TABLE>" >> $TSF_TEST/tsftest.html
echo "" >>  $TSF_TEST/tsftest.html
}

write_html_footer() {
echo " </BODY>" >> $TSF_TEST/tsftest.html
echo "</HTML>" >> $TSF_TEST/tsftest.html
return 0;
}

write_entry() {
echo "     <TR> " >> $TSF_TEST/tsftest.html
echo "       <TD>$1</TD>" >> $TSF_TEST/tsftest.html 
case "$2" in
   "PASSED") echo "       <TD><FONT COLOR=#00FF00>$2</FONT></TD> " >> $TSF_TEST/tsftest.html;;
   "FAILED") echo "       <TD><FONT COLOR=#FF0000>$2</FONT></TD> " >> $TSF_TEST/tsftest.html;;
   "EXCEPTION") echo "       <TD><FONT COLOR=#FF0000>$2</FONT></TD> " >> $TSF_TEST/tsftest.html;;
   "CRASHED") echo "       <TD><FONT COLOR=#FF0000>CRASHED</TD> " >> $TSF_TEST/tsftest.html;;
          *) echo "       <TD><FONT COLOR=#FF0000>UNKNOWN</FONT></TD> " >> $TSF_TEST/tsftest.html;;

    esac
echo "       <TD>$3</TD> " >> $TSF_TEST/tsftest.html
echo "     </TR> " >> $TSF_TEST/tsftest.html

return 0;
}
#################################################################################
#
# Set defaults and clean up any old files
#
#################################################################################
TSF_ROOT=$TRILINOS_HOME/packages/TSF
TSF_TEST=$TSF_ROOT/examples/Tests

# /bin/rm -rf cvsusers

#################################################################################
# Parse options
#################################################################################
#while getopts co: OPTION
#do
#case $OPTION in
#    c) CHECKOUT=1;;
#    o) OUTPUT=$OPTARG;;
#esac
#done
#
# if ([ $CHECKOUT -eq 1 ])
# then
#################################################################################
#
# Check out TSF, and make tests
#
#################################################################################
#    export CVS_RSH=ssh
#    cvsdir=":ext:mwboldt@iterative.ca.sandia.gov:/usr/local/CVS"
#    /bin/rm -rf ./$TSF_ROOT
#    cvs -d $cvsdir co $TSF_ROOT
#    (cd $TSF_ROOT; ./configure; make tests >& maketests.log)
# fi
#

# for file in fitSine cgTest findMin petraPoissonTest powerMethodTest serialSolverTest
# do
#     ./$file.e &> $file.log
#     cat $file.log >> $TSF_TEST/tsftest.out
#     echo -n $file >> $TSF_TEST/results.csv
#     if [ -s $file.csv ]
#     then
# 	cat $file.csv >> $TSF_TEST/results.csv
#     else
# 	echo " CRASHED see_log" >> $TSF_TEST/results.csv
#     fi
# done
# cd $TSF_TEST
#################################################################################
# Extract all cvs users who executed a commit yesterday
# Create html file with results
#################################################################################
#cvs history -a -c -D "1 day ago" | cut -c26-32 | uniq > cvsusers

if [ -e tsftest.html ]
then
    mv tsftest.html old.html
fi
write_html_header
write_table_header `date -Iminutes`

while read testname teststatus testcomments
do
    write_entry $testname $teststatus $testcomments
done <results.csv
write_table_footer
if [ -e old.html ]
then
    # Combine new and old results leaving only one html header
    tail -n +7 old.html >> tsftest.html
    rm -f old.html
else
    # Only write html footer if these are the first results
    write_html_footer
fi

exit 0
