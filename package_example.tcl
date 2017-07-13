# Create the namespace
namespace eval ::HelloWorld {
 
	# Export proc
    # namespace export $::HelloWorld
    namespace export talk
 
	# Set up Variable
	set version 1.0
	set HelloWorldDescription "HelloWorld"
 
    # Variable for the path of the script
    variable home [file join [pwd] [file dirname [info script]]]
}
 
# Definition of the procedure talk
proc ::HelloWorld::TestRun {} {
	puts $HelloWorld::HelloWorldDescription
}
# Definition of the procedure getpath
proc ::HelloWorld::GetPath {} {
	variable home
        return $home
}

proc talk {args} \
{
	return [eval ::HelloWorld::TestRun $args]
}
package provide HelloWorld $HelloWorld::version
# package require Tcl 8.4