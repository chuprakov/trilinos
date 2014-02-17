# awk script to provide context for MueLu timers
#
# Don't run this script directly, instead use compareTE.sh.
#
# Note: This requires GNU awk and won't run in compatibility mode.
#
# Note: You might need to define AWKPATH if mueprof.sh can't find this script.
#
# Thanks to http://www.grymoire.com/Unix/Awk.html and CMS for help.

###############################################################################
#stuff that happens before processing any files
BEGIN {
    whichBlockToAnalyze=blockNumber; #passed in from shell script
    foundTimerBlock=0;
    foundTimersToReport=0;
    possibleTotalLabels[1] = "ScalingTest: 2 - MueLu Setup";
    possibleTotalLabels[2] = "MueLu: Hierarchy: Setup [(]total[)]";
    possibleTotalLabels[3] = "nalu MueLu preconditioner setup";
    tt=""
}

###############################################################################
#stuff that happens when processing all files
{
    if (tt != FILENAME) {
      tt=FILENAME
      #reset for next file
      whichBlockToAnalyze=blockNumber;
      foundTimerBlock=0;
      foundTimersToReport=0;
      #hack, since we can't count on having "Linear algebra library: Tpetra" print in nalu
      linalg[FILENAME] = "Epetra"
    }

    #indicates start of MueLu timer output
    if (match($0,"^Timer Name")) {
      #logic to ensure only the ith solver block is analyzed
      whichBlockToAnalyze--;
      if (whichBlockToAnalyze>0)
        foundTimerBlock=0;
      if (whichBlockToAnalyze==0)
        foundTimerBlock=1;
      if (whichBlockToAnalyze<0) {
        nextfile; #we've moved past the solver block of interest,
                  #stop analyzing this file
      }
    }

    if (foundTimerBlock) {
      if (match($0,"^MueLu: ")) {
        # matched a MueLu timer
        #fix minor difference between Epetra and Tpetra smoother tags
        sub(/Ifpack2Smoother/,"IfpackSmoother");
        sub(/Amesos2Smoother/,"AmesosSmoother");
        if (match($0,"[(]level=[0-9][)]")) {
          # timer is level-specific (and by its nature excludes calls to child factories)
          foundTimersToReport=1;
          factAndLevel = substr($0,1,RSTART-1+RLENGTH);
          alltimes = substr($0,RSTART+RLENGTH);
          cutCmd="cut -f3 -d')' | cut -f1 -d'('"
          maxtime = ExtractTime(alltimes,cutCmd);
          if (match(factAndLevel,"MueLu: Hierarchy: Solve")) {
            #TODO figure out which solve labels to pull out
            solveLabels[factAndLevel] = factAndLevel;
            solveTimes[factAndLevel,FILENAME] = maxtime;
          } else {
            setupLabels[factAndLevel] = factAndLevel;
            setupTimes[factAndLevel,FILENAME] = maxtime;
          }
        }
      }

      if (match($0,"^YY Tpetra Transpose Only")) {
        linalg[FILENAME] = "Tpetra"
      }

      # Check for any reported total setup time.  This is printed as a sanity check
      # against the running total.
      for (i in possibleTotalLabels) {
        if (match($0,possibleTotalLabels[i])) {
          pattern = substr($0,RSTART,RLENGTH);
          alltimes = substr($0,RSTART+RLENGTH);
          cutCmd="cut -f3 -d')' | cut -f1 -d'('"
          TotalSetup[pattern,FILENAME] = ExtractTime(alltimes,cutCmd);
        }
      }
    } #if (foundTimerBlock)

}

###############################################################################
function PrintHeader(description,linalg)
{
  space = " ";
  printf("%60s      ",toupper(description));
  for (j in linalg) {
    printf("%10s  ",linalg[j]);
    printf(" (total)");
  }
  if (length(linalg)>1)
    printf("%13s  ","T/E ratio");
  printf("\n%60s          -------------------------------------------------\n",space);
}

###############################################################################
# For sanity purposes, print the total that MueLu prints
function PrintTotalSetupTime()
{
  printf("%60s          ----------------\n"," ");
  for (i in TotalSetup) {
    split(i,sep,SUBSEP); #breaks multiarray index i up into its constituent parts.
                         #we only want the first one, sep[1]
    printf("%60s  ==>   %6.3f seconds       \n",sep[1],TotalSetup[i]);
  }
}

###############################################################################
function SortAndReportTimings(libToSortOn, timerLabels, timerValues, linalg)
{
  if (length(linalg) == 1)
    for (i in linalg) sortByLib = linalg[i]
  len=0
  for (i in timerLabels) {
    valuesAndLabels[len] = timerLabels[i];
    for (j in linalg)
      valuesAndLabels[len] = valuesAndLabels[len] "@" timerValues[i,j]
    len++
  }

  #combine timer label and corresponding values into one array entry
  sortField=2
  for (i in linalg) {
    if (linalg[i] == libToSortOn) break
    sortField++
  }
  sortCmd = "sort -k" sortField" -n -t @"
  SystemSort(valuesAndLabels,sortedValuesAndLabels,sortCmd);

  jj=2
  for (i in linalg) {
    runningTotals[i] = 0
    if (linalg[i] == "Epetra") epetraInd = jj;
    if (linalg[i] == "Tpetra") tpetraInd = jj;
    jj++
  }
  for (j=1; j<=len; j++) {
    split(sortedValuesAndLabels[j],fields,"@");
    printf("%60s  ==> ",fields[1]);
    offset=2
    for (i in linalg) {
      runningTotals[i] += fields[offset]
      printf("%10.2f   (%5.1f)",fields[offset],runningTotals[i]);
      offset++
    }
    if (fields[epetraInd] != 0 && length(linalg) > 1) {
      printf("%11.2f   ",fields[tpetraInd] / fields[epetraInd]);
    }
    printf("\n");
  }

}
###############################################################################

function SystemSort(arrayToSort, sortedArray, sortCmd)
{
  # Write to the coprocess...
  for(i in arrayToSort) print arrayToSort[i] |& sortCmd
  close(sortCmd, "to")
  i=0
  # ...and read back from it
  while((sortCmd |& getline sortedArray[++i]) > 0)
    ;
  close(sortCmd)
}

###############################################################################

function ExtractTime(alltimes,cutCmd)
{
  alltimes = RemoveWhiteSpace(alltimes);
  # system call to cutCmd
  # see GNU Awk manual section 12.3, Two-Way Communications with Another Process
  print alltimes |& cutCmd
  close(cutCmd,"to")
  cutCmd |& getline themax
  close(cutCmd)
  themax = RemoveWhiteSpace(themax);
  return (themax);
}

###############################################################################

function RemoveWhiteSpace(theString)
{
  sub(/^[ ]*/,"",theString); sub(/[ ]*$/,"",theString);
  return (theString)
}
###############################################################################

###############################################################################
#stuff that happens after processing all files
###############################################################################
END {

  if (foundTimersToReport) {
    PrintHeader("Setup times (level specific) excluding child calls ",linalg);
    SortAndReportTimings(sortByLib,setupLabels,setupTimes,linalg);
    PrintTotalSetupTime();
  } else {
    printf("\nFound only %d solver blocks, you requested block %d.\n",blockNumber-whichBlockToAnalyze,blockNumber);
  }
  printf("\n");

}
