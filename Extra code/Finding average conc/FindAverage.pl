#!/usr/bin/perl

# uses data point locations to place and scale an initial mesh

use strict;

my ($node,$concen);


################################################
############## INPUTS #################

my $nodefile1  = 'HealthyTimeNoOutletResults_0.002.part0.exnode';
my $nodefile2  = 'HealthyTimeNoOutletResults_0.004.part0.exnode';
my $nodefile3  = 'HealthyTimeNoOutletResults_0.006.part0.exnode';
my $nodefile4  = 'HealthyTimeNoOutletResults_0.008.part0.exnode';
my $nodefile5  = 'HealthyTimeNoOutletResults_0.01.part0.exnode';
my $nodefile6  = 'HealthyTimeNoOutletResults_0.012.part0.exnode';
#my $nodefile7  = 'HealthyTimeResults_0.7.part0.exnode';
#my $nodefile8  = 'HealthySteadyStateResults.part0.exnode';
my $nodefile7  = 'HealthySteadyStateResults.part0.exnode';

my $time1 = 0.002;
my $time2 = 0.004;
my $time3 = 0.006;
my $time4 = 0.008;
my $time5 = 0.01;
my $time6 = 0.012;
#my $time7 = 0.7;
#my $time8 = 0.8;
my $time7 = 1;

#######################################################


print "open initial node file $nodefile1 \n";

open EXNODE, "<$nodefile1" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;

my $sumconc =0.0;

my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Node:\s+(\d+)/) {
	$node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++) { ## for x,y,z
	    $line = <EXNODE>;
        }
	    $line = <EXNODE>;  

            $concen = $line;
            ## print "t1.  $concen  \n";
            $sumconc = $sumconc + $concen;
	}
	
    $line = <EXNODE>;    
}

close EXNODE;

my $avg = $sumconc/$NumberofNodes;

print "number of nodes =   $NumberofNodes \n";
print "\n";
print "------printing average-------\n";
print "$time1  $avg  \n";

#######################################################################



print "\n";
print "open initial node file $nodefile2 \n";

open EXNODE, "<$nodefile2" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;

$sumconc =0.0;


my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Node:\s+(\d+)/) {
	$node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++) { ## for x,y,z
	    $line = <EXNODE>;
        }
	    $line = <EXNODE>;  

            $concen = $line;
            ## print "t1.  $concen  \n";
            $sumconc = $sumconc + $concen;
	}
	
    $line = <EXNODE>;    
}

close EXNODE;

my $avg = $sumconc/$NumberofNodes;

print "\n";
print "------printing average-------\n";
print "$time2  $avg  \n";


#############################################################################



print "\n";
print "open initial node file $nodefile3 \n";

open EXNODE, "<$nodefile3" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;

$sumconc =0.0;


my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Node:\s+(\d+)/) {
	$node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++) { ## for x,y,z
	    $line = <EXNODE>;
        }
	    $line = <EXNODE>;  

            $concen = $line;
            ## print "t1.  $concen  \n";
            $sumconc = $sumconc + $concen;
	}
	
    $line = <EXNODE>;    
}

close EXNODE;

my $avg = $sumconc/$NumberofNodes;

print "\n";
print "------printing average-------\n";
print "$time3  $avg  \n";


#############################################################################



print "\n";
print "open initial node file $nodefile4 \n";

open EXNODE, "<$nodefile4" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;

$sumconc =0.0;


my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Node:\s+(\d+)/) {
	$node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++) { ## for x,y,z
	    $line = <EXNODE>;
        }
	    $line = <EXNODE>;  

            $concen = $line;
            ## print "t1.  $concen  \n";
            $sumconc = $sumconc + $concen;
	}
	
    $line = <EXNODE>;    
}

close EXNODE;

my $avg = $sumconc/$NumberofNodes;

print "\n";
print "------printing average-------\n";
print "$time4  $avg  \n";


#############################################################################



print "\n";
print "open initial node file $nodefile5 \n";

open EXNODE, "<$nodefile5" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;

$sumconc =0.0;


my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Node:\s+(\d+)/) {
	$node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++) { ## for x,y,z
	    $line = <EXNODE>;
        }
	    $line = <EXNODE>;  

            $concen = $line;
            ## print "t1.  $concen  \n";
            $sumconc = $sumconc + $concen;
	}
	
    $line = <EXNODE>;    
}

close EXNODE;

my $avg = $sumconc/$NumberofNodes;

print "\n";
print "------printing average-------\n";
print "$time5  $avg  \n";


#############################################################################



print "\n";
print "open initial node file $nodefile6 \n";

open EXNODE, "<$nodefile6" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;

$sumconc =0.0;


my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Node:\s+(\d+)/) {
	$node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++) { ## for x,y,z
	    $line = <EXNODE>;
        }
	    $line = <EXNODE>;  

            $concen = $line;
            ## print "t1.  $concen  \n";
            $sumconc = $sumconc + $concen;
	}
	
    $line = <EXNODE>;    
}

close EXNODE;

my $avg = $sumconc/$NumberofNodes;

print "\n";
print "------printing average-------\n";
print "$time6  $avg  \n";


#############################################################################



print "\n";
print "open initial node file $nodefile7 \n";

open EXNODE, "<$nodefile7" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;

$sumconc =0.0;


my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Node:\s+(\d+)/) {
	$node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++) { ## for x,y,z
	    $line = <EXNODE>;
        }
	    $line = <EXNODE>;  

            $concen = $line;
            ## print "t1.  $concen  \n";
            $sumconc = $sumconc + $concen;
	}
	
    $line = <EXNODE>;    
}

close EXNODE;

my $avg = $sumconc/$NumberofNodes;

print "\n";
print "------printing average-------\n";
print "$time7  $avg  \n";


#############################################################################



#print "\n";
#print "open initial node file $nodefile8 \n";

#open EXNODE, "<$nodefile8" or die "\033[31mError: Can't open node file\033[0m ";
#my $NumberofNodes=0;

#$sumconc =0.0;


#my $line = <EXNODE>;
#while ($line) {

#    if ($line =~ /Node:\s+(\d+)/) {
#	$node = $1;
#	$NumberofNodes=$NumberofNodes+1;
#    
#	for (my $i=1;$i<=3;$i++) { ## for x,y,z
#	    $line = <EXNODE>;
#        }
#	    $line = <EXNODE>;  

#            $concen = $line;
#            ## print "t1.  $concen  \n";
#            $sumconc = $sumconc + $concen;
#	}
#	
#    $line = <EXNODE>;    
#}

#close EXNODE;

#my $avg = $sumconc/$NumberofNodes;

#print "\n";
#print "------printing average-------\n";
#print "$time8  $avg  \n";


##############################################################################




