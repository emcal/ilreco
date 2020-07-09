#!/usr/bin/tclsh

set fidc [open oc]
set fidf [open of]

while {([eof $fidc] == 0) && ([eof $fidf] == 0)} {
  set linec [gets $fidc]
  set linef [gets $fidf]
  if {([llength $linec] != [llength $linef])} {
    puts $linec
    puts $linef
    puts "lines mimatch"
    break
  }
  set dif 0
  for {set i 0} {($i < [llength $linec])} {incr i} {
    set vf [lindex $linef $i]
    set vc [lindex $linec $i]
    if [catch {set dummy [expr $vc+$vf]}] {
      if {([catch {set dummy [expr $vc]}] == 0) && \
          ([catch {set dummy [expr $vf]}] == 0)} {set dif 1}
      break
    }
    if {($vc != $vf)} {
#     if {([expr abs($vc - $vf)]) > 3.3e-5} {set dif 1}
      if {([expr abs($vc+$vf)] != 0.)} {
        if {([expr abs($vc-$vf)/abs($vc+$vf)] > 0.05)} {set dif 1}
      }
    }
  }
  if $dif {
    puts $linec
    puts $linef
    puts ""
  }
}

close $fidc
close $fidf
