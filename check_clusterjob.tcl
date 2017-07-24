# CheckCluster job
# check SGE array task on TCBG local cluster
# Author: Maximilian Scheurer, April 2016, Mail: mscheurer@ks.uiuc.edu
namespace eval ::CheckCluster {
	set version 0.1
	set description "check array job status on TCBG cluster"
}
package provide CheckCluster $CheckCluster::version

proc check_clusterjob {args} \
{
	return [eval ::CheckCluster::check_clusterjob $args]
}

proc ::CheckCluster::check_clusterjob {uname jobid ntasks} \
{
	set output [exec qstat -u $uname]
	set checkjob_rc [catch { exec qstat -j $jobid } msg]
	if {$checkjob_rc == 1} {
		puts "job is not in the queue"		
		return 0
	} else {
		if {$msg != ""} {
			set currentout [open "current_stat.txt" w]
			puts $currentout $output
			close $currentout
			set column_extract [exec awk "/$jobid/ {print \$1 \"\t\" \$5 \"\t\" \$(NF)}" current_stat.txt]
			set lines [split $column_extract "\n"]
			set n_queue 0
			set njobs $ntasks
			set running_jobs []
			foreach line $lines {
				set status [lindex $line 1]
				switch -exact -- $status {
					"r" {
						lappend running_jobs [lindex $line 2]
					}	
					"t" {
						lappend running_jobs [lindex $line 2]
					}
					"qw" {
						set inf [lindex $line 2]
						set parse [split $inf ":"]
						set range [lindex $parse 0]
						if {[string  first "," $range] != -1} {
							set idx1 [lindex [split $range ","] 0]
							set idx2 [lindex [split $range ","] 1]
						} else {
							set idx1 [lindex [split $range "-"] 0]
							set idx2 [lindex [split $range "-"] 1]
						}
						set n_queue [expr $idx2 - $idx1 + 1]
					}
					"hqw" {
						set inf [lindex $line 2]
						set parse [split $inf ":"]
						set range [lindex $parse 0]
						if {[string  first "," $range] != -1} {
							set idx1 [lindex [split $range ","] 0]
							set idx2 [lindex [split $range ","] 1]
						} else {
							set idx1 [lindex [split $range "-"] 0]
							set idx2 [lindex [split $range "-"] 1]
						}
						set n_queue [expr $idx2 - $idx1 + 1]
					}
					default {
				
					}
				}
			}
			set n_running [expr [lindex $running_jobs end]-[lindex $running_jobs 0] + 1]
			set pc_queued [expr double($n_queue)/double($njobs) * 100]
			set pc_running [expr double($n_running)/double($njobs) * 100]
			set n_fin [expr $njobs - $n_running - $n_queue]
			set pc_fin [expr double($n_fin)/double($njobs) * 100]
			puts "running: $n_running ($pc_running\%), queued: $n_queue ($pc_queued\%), finished: $n_fin ($pc_fin\%)"
			return 1
		} else {
			return 0
		}
	}
}
