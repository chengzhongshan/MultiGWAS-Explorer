#!/usr/bin/perl -w
use strict;
use FileHandle;
my $exclude_jobs_ref;
if (@ARGV==0) {
   $exclude_jobs_ref=[0];
}else{
   $exclude_jobs_ref=\@ARGV;  
}
#Noeed to exclude the PID for the current analysis of ps and grep;
my $FH=FileHandle->new("ps a |grep 'usr\/bin\/(ps|grep)' -Pv|");#ps ax or other commands;
my @jobs;
while (my $pid=<$FH>) {
   chomp($pid);
   $pid=~s/^\s+//;
   $pid=~s/ +/\t/g;
   if ($pid=~/^\d+/) {
       my @as=split("\t| ",$pid);
       my $job=$as[0];
       #print $job,"\n";
       my $t=0;
       foreach my $j(@$exclude_jobs_ref){
          $t++ if ($job==$j);
       }
       push @jobs,$job if $t==0;
       #print $pid,"\n";
   }
}
$FH->close;
print join(" ",@jobs),"\n" if @jobs>0;




