# ======================================================================
# Cylinder axial projection to silicone outer surfaces
# HyperMesh 2025
#
# Workflow:
#   1. Select silicone outer surfaces once.
#   2. Select one cylinder Solid.
#   3. Manually define the cylinder axis.
#   4. Search surface intersections along the cylinder axis.
#   5. Create a node at the projected surface point.
#   6. Repeat for the next cylinder.
#   7. Export all results to CSV.
#
# System ID:
#   System ID = 1 is used for local-coordinate output.
#
# Projection modes:
#   closest  : Search both axis directions and use the nearest point.
#   positive : Only use the positive direction defined in the direction panel.
#   negative : Only use the opposite direction.
#
# Result status:
#   EXACT            : Exact line/surface intersection.
#   APPROX_NEAR_MISS : No exact intersection, but the miss distance was
#                      smaller than near_miss_tolerance.
#   NO_HIT           : No acceptable intersection.
# ======================================================================


namespace eval CYLAX {

    variable cfg

    # --------------------------------------------------------------
    # User settings
    # --------------------------------------------------------------

    array set cfg [list \
        system_id                  1 \
        projection_mode           "closest" \
        create_projection_node    1 \
        export_csv                1 \
        prompt_after_each         1 \
        keep_failed_axis_line     1 \
        allow_near_miss_fallback  1 \
        near_miss_tolerance       0.01 \
        minimum_projection_dist   1.0e-7 \
        ray_length_scale          1.50 \
        manual_ray_half_length    0.0 \
        decimal_places            6 \
    ]

    variable results {}
}


# ======================================================================
# Vector functions
# ======================================================================

proc CYLAX::VSub {a b} {

    return [list \
        [expr {[lindex $a 0] - [lindex $b 0]}] \
        [expr {[lindex $a 1] - [lindex $b 1]}] \
        [expr {[lindex $a 2] - [lindex $b 2]}]]
}


proc CYLAX::VAddScaled {point direction scale} {

    return [list \
        [expr {[lindex $point 0] +
               $scale*[lindex $direction 0]}] \
        [expr {[lindex $point 1] +
               $scale*[lindex $direction 1]}] \
        [expr {[lindex $point 2] +
               $scale*[lindex $direction 2]}]]
}


proc CYLAX::VDot {a b} {

    return [expr {
        [lindex $a 0]*[lindex $b 0] +
        [lindex $a 1]*[lindex $b 1] +
        [lindex $a 2]*[lindex $b 2]
    }]
}


proc CYLAX::VNorm {a} {

    return [expr {sqrt([CYLAX::VDot $a $a])}]
}


proc CYLAX::VUnit {a} {

    set length [CYLAX::VNorm $a]

    if {$length < 1.0e-15} {
        error "The direction vector has zero length."
    }

    return [list \
        [expr {[lindex $a 0]/$length}] \
        [expr {[lindex $a 1]/$length}] \
        [expr {[lindex $a 2]/$length}]]
}


# ======================================================================
# Formatting and coordinate conversion
# ======================================================================

proc CYLAX::Fmt {value digits} {

    return [format "%.${digits}f" $value]
}


proc CYLAX::FmtPoint {point digits} {

    return [list \
        [CYLAX::Fmt [lindex $point 0] $digits] \
        [CYLAX::Fmt [lindex $point 1] $digits] \
        [CYLAX::Fmt [lindex $point 2] $digits]]
}


proc CYLAX::ToLocal {system_id point} {

    set x [lindex $point 0]
    set y [lindex $point 1]
    set z [lindex $point 2]

    return [list \
        [hm_xpointlocal $system_id $x $y $z] \
        [hm_ypointlocal $system_id $x $y $z] \
        [hm_zpointlocal $system_id $x $y $z]]
}


# ======================================================================
# Direction-panel result handling
# ======================================================================

proc CYLAX::ParseDirection {raw_direction} {

    if {[llength $raw_direction] == 0} {
        error "No direction was defined."
    }

    # Compatibility with nested return formats.
    if {
        [llength $raw_direction] == 1 &&
        [llength [lindex $raw_direction 0]] >= 3
    } {
        set raw_direction [lindex $raw_direction 0]
    }

    if {[llength $raw_direction] < 3} {
        error "Invalid direction data: $raw_direction"
    }

    set direction [list \
        [lindex $raw_direction 0] \
        [lindex $raw_direction 1] \
        [lindex $raw_direction 2]]

    return [CYLAX::VUnit $direction]
}


# ======================================================================
# Direction filtering
# ======================================================================

proc CYLAX::DirectionAllowed {
    signed_distance projection_mode minimum_distance
} {

    switch -- $projection_mode {

        "closest" {
            return [expr {
                abs($signed_distance) > $minimum_distance
            }]
        }

        "positive" {
            return [expr {
                $signed_distance > $minimum_distance
            }]
        }

        "negative" {
            return [expr {
                $signed_distance < -$minimum_distance
            }]
        }

        default {
            error "Unknown projection_mode: $projection_mode"
        }
    }
}


proc CYLAX::CandidateKey {signed_distance projection_mode} {

    switch -- $projection_mode {

        "closest" {
            return [expr {abs($signed_distance)}]
        }

        "positive" {
            return $signed_distance
        }

        "negative" {
            return [expr {-$signed_distance}]
        }

        default {
            error "Unknown projection_mode: $projection_mode"
        }
    }
}


# ======================================================================
# Axis-line management
# ======================================================================

proc CYLAX::DeleteLine {line_id} {

    if {$line_id eq "" || $line_id <= 0} {
        return
    }

    catch {
        *createmark lines 1 $line_id
        *deletemark lines 1
    }
}


proc CYLAX::CreateAxisLine {
    centroid direction half_length
} {

    set start_point \
        [CYLAX::VAddScaled $centroid $direction \
            [expr {-$half_length}]]

    set end_point \
        [CYLAX::VAddScaled $centroid $direction \
            $half_length]

    set x1 [lindex $start_point 0]
    set y1 [lindex $start_point 1]
    set z1 [lindex $start_point 2]

    set x2 [lindex $end_point 0]
    set y2 [lindex $end_point 1]
    set z2 [lindex $end_point 2]

    *linecreatestraight $x1 $y1 $z1 $x2 $y2 $z2

    set line_id [hm_latestentityid lines]

    if {$line_id eq "" || $line_id <= 0} {
        error "Failed to create the temporary cylinder-axis line."
    }

    return $line_id
}


# ======================================================================
# Automatic line length
# ======================================================================

proc CYLAX::GetAxisHalfLength {
    centroid bbox scale manual_length
} {

    if {$manual_length > 0.0} {
        return $manual_length
    }

    foreach {
        xmin ymin zmin
        xmax ymax zmax
    } [lrange $bbox 0 5] {}

    set bbox_center [list \
        [expr {($xmin+$xmax)/2.0}] \
        [expr {($ymin+$ymax)/2.0}] \
        [expr {($zmin+$zmax)/2.0}]]

    set half_diagonal [expr {
        0.5*sqrt(
            ($xmax-$xmin)*($xmax-$xmin) +
            ($ymax-$ymin)*($ymax-$ymin) +
            ($zmax-$zmin)*($zmax-$zmin)
        )
    }]

    set centroid_to_box_center \
        [CYLAX::VNorm \
            [CYLAX::VSub $centroid $bbox_center]]

    set half_length [expr {
        $scale*(
            $centroid_to_box_center +
            $half_diagonal
        )
    }]

    if {$half_length <= 0.0} {
        set half_length 1.0
    }

    return $half_length
}


# ======================================================================
# Find the cylinder-axis intersection with the selected surfaces
#
# Return dictionary fields:
#   status
#   line_id
#   surface_id
#   point
#   signed_distance
#   miss_distance
#   nearest_miss_surface
#   nearest_miss_distance
# ======================================================================

proc CYLAX::FindProjection {
    centroid
    direction
    half_length
    surface_ids
    projection_mode
    minimum_distance
    allow_near_miss
    near_miss_tolerance
    keep_failed_line
} {

    set line_id \
        [CYLAX::CreateAxisLine \
            $centroid $direction $half_length]

    set best_exact {}
    set best_exact_key 1.0e300

    set best_miss {}
    set best_miss_distance 1.0e300
    set best_miss_surface -1

    foreach surface_id $surface_ids {

        if {[catch {
            set data \
                [hm_getclosestpointsbetweenlinesurface \
                    $line_id $surface_id]
        } query_error]} {

            puts "  Surface $surface_id query error: $query_error"
            continue
        }

        set value_count [llength $data]

        # ----------------------------------------------------------
        # No exact intersection:
        #
        # xs ys zs = closest surface point
        # xc yc zc = closest line point
        # dist      = minimum distance
        # ----------------------------------------------------------

        if {$value_count == 7} {

            foreach {
                xs ys zs
                xc yc zc
                miss_distance
            } $data {}

            set line_point \
                [list $xc $yc $zc]

            set signed_distance \
                [CYLAX::VDot \
                    [CYLAX::VSub $line_point $centroid] \
                    $direction]

            if {
                [CYLAX::DirectionAllowed \
                    $signed_distance \
                    $projection_mode \
                    $minimum_distance]
            } {

                if {$miss_distance < $best_miss_distance} {

                    set best_miss_distance $miss_distance
                    set best_miss_surface $surface_id

                    set best_miss [dict create \
                        status          "APPROX_NEAR_MISS" \
                        line_id         $line_id \
                        surface_id      $surface_id \
                        point           [list $xs $ys $zs] \
                        signed_distance $signed_distance \
                        miss_distance   $miss_distance]
                }
            }

            continue
        }


        # ----------------------------------------------------------
        # Exact intersection:
        # Each intersection contains six values:
        #
        # xs ys zs xc yc zc
        # ----------------------------------------------------------

        if {
            $value_count < 6 ||
            [expr {$value_count % 6}] != 0
        } {
            puts "  Surface $surface_id returned unexpected data:"
            puts "  $data"
            continue
        }

        for {set i 0} {$i < $value_count} {incr i 6} {

            set xs [lindex $data $i]
            set ys [lindex $data [expr {$i+1}]]
            set zs [lindex $data [expr {$i+2}]]

            set xc [lindex $data [expr {$i+3}]]
            set yc [lindex $data [expr {$i+4}]]
            set zc [lindex $data [expr {$i+5}]]

            set line_point [list $xc $yc $zc]

            set signed_distance \
                [CYLAX::VDot \
                    [CYLAX::VSub $line_point $centroid] \
                    $direction]

            if {
                ![CYLAX::DirectionAllowed \
                    $signed_distance \
                    $projection_mode \
                    $minimum_distance]
            } {
                continue
            }

            set candidate_key \
                [CYLAX::CandidateKey \
                    $signed_distance \
                    $projection_mode]

            if {$candidate_key < $best_exact_key} {

                set best_exact_key $candidate_key

                set best_exact [dict create \
                    status          "EXACT" \
                    line_id         $line_id \
                    surface_id      $surface_id \
                    point           [list $xs $ys $zs] \
                    signed_distance $signed_distance \
                    miss_distance   0.0]
            }
        }
    }


    # --------------------------------------------------------------
    # Exact result has priority
    # --------------------------------------------------------------

    if {[dict size $best_exact] > 0} {

        CYLAX::DeleteLine $line_id

        dict set best_exact line_id ""

        return $best_exact
    }


    # --------------------------------------------------------------
    # Optional small near-miss fallback
    # --------------------------------------------------------------

    if {
        $allow_near_miss &&
        [dict size $best_miss] > 0 &&
        $best_miss_distance <= $near_miss_tolerance
    } {

        CYLAX::DeleteLine $line_id

        dict set best_miss line_id ""

        return $best_miss
    }


    # --------------------------------------------------------------
    # No acceptable point
    # --------------------------------------------------------------

    if {!$keep_failed_line} {
        CYLAX::DeleteLine $line_id
        set line_id ""
    }

    return [dict create \
        status                "NO_HIT" \
        line_id               $line_id \
        surface_id            -1 \
        point                 {} \
        signed_distance       "" \
        miss_distance         "" \
        nearest_miss_surface  $best_miss_surface \
        nearest_miss_distance $best_miss_distance]
}


# ======================================================================
# Continue dialog
# ======================================================================

proc CYLAX::AskContinue {solid_id} {

    variable cfg

    if {!$cfg(prompt_after_each)} {
        return 1
    }

    if {[catch {
        set answer [tk_messageBox \
            -title "Cylinder axial projection" \
            -message \
                "Solid $solid_id has been processed.\n\nContinue with the next cylinder?" \
            -type yesno \
            -icon question]
    }]} {
        return 1
    }

    return [expr {$answer eq "yes"}]
}


# ======================================================================
# CSV output
# ======================================================================

proc CYLAX::ExportCSV {} {

    variable cfg
    variable results

    if {!$cfg(export_csv)} {
        return
    }

    if {[llength $results] <= 1} {
        puts "No projection results are available for CSV export."
        return
    }

    set csv_path [tk_getSaveFile \
        -title "Save cylinder axial projection results" \
        -defaultextension ".csv" \
        -filetypes {
            {"CSV files" {.csv}}
            {"All files" {*}}
        } \
        -initialfile \
            "Cylinder_Axial_Projection_System_1.csv"]

    if {$csv_path eq ""} {
        puts "CSV export was cancelled."
        return
    }

    set file_handle ""

    if {[catch {

        set file_handle [open $csv_path "w"]

        fconfigure $file_handle \
            -encoding utf-8 \
            -translation lf

        # UTF-8 BOM for Excel.
        puts -nonewline $file_handle "\ufeff"

        foreach row $results {
            puts $file_handle [join $row ","]
        }

        close $file_handle
        set file_handle ""

    } write_error]} {

        if {$file_handle ne ""} {
            catch {close $file_handle}
        }

        puts "CSV export failed:"
        puts $write_error
        return
    }

    puts ""
    puts "CSV saved to:"
    puts $csv_path
}


# ======================================================================
# Main procedure
# ======================================================================

proc CYLAX::Main {} {

    variable cfg
    variable results

    puts ""
    puts "============================================================"
    puts "Cylinder axial projection - HyperMesh 2025"
    puts "System ID              = $cfg(system_id)"
    puts "Projection mode        = $cfg(projection_mode)"
    puts "Near-miss fallback     = $cfg(allow_near_miss_fallback)"
    puts "Near-miss tolerance    = $cfg(near_miss_tolerance)"
    puts "Keep failed axis lines = $cfg(keep_failed_axis_line)"
    puts "============================================================"


    # --------------------------------------------------------------
    # Validate local system
    # --------------------------------------------------------------

    if {[catch {
        hm_xpointlocal $cfg(system_id) 0.0 0.0 0.0
    } system_error]} {

        puts "ERROR: Cannot use System ID $cfg(system_id)."
        puts $system_error
        return
    }


    # --------------------------------------------------------------
    # Select silicone outer surfaces once
    # --------------------------------------------------------------

    *createmarkpanel surfs 2 \
        "Select all silicone OUTER surfaces, then confirm"

    set surface_ids [hm_getmark surfs 2]

    if {[llength $surface_ids] == 0} {
        puts "No silicone surfaces were selected."
        return
    }

    puts "Selected silicone surfaces: [llength $surface_ids]"


    # --------------------------------------------------------------
    # Get silicone bounding box
    # --------------------------------------------------------------

    if {[catch {
        set bbox [hm_getboundingbox surfs 2 0 0 0]
    } bbox_error]} {

        puts "ERROR: Cannot calculate silicone bounding box."
        puts $bbox_error
        return
    }

    if {[llength $bbox] < 6} {
        puts "ERROR: Invalid silicone bounding box."
        return
    }


    # --------------------------------------------------------------
    # Results header
    # --------------------------------------------------------------

    set results [list [list \
        "Solid_ID" \
        "Status" \
        "Centroid_Global_X" \
        "Centroid_Global_Y" \
        "Centroid_Global_Z" \
        "Centroid_Local_X" \
        "Centroid_Local_Y" \
        "Centroid_Local_Z" \
        "Axis_Global_X" \
        "Axis_Global_Y" \
        "Axis_Global_Z" \
        "Projection_Side" \
        "Surface_ID" \
        "Projected_Global_X" \
        "Projected_Global_Y" \
        "Projected_Global_Z" \
        "Projected_Local_X" \
        "Projected_Local_Y" \
        "Projected_Local_Z" \
        "Signed_Axis_Distance" \
        "Absolute_Axis_Distance" \
        "Miss_Distance" \
        "Created_Node_ID" \
        "Debug_Line_ID"]]


    set success_count 0
    set approximate_count 0
    set failed_count 0


    # --------------------------------------------------------------
    # Process cylinders one by one
    # --------------------------------------------------------------

    while {1} {

        puts ""
        puts "Select one cylinder Solid."
        puts "Cancel or confirm an empty selection to finish."

        *createmarkpanel solids 1 \
            "Select ONE cylinder Solid; empty selection finishes"

        set solid_ids [hm_getmark solids 1]

        if {[llength $solid_ids] == 0} {
            break
        }

        if {[llength $solid_ids] != 1} {
            puts "Please select exactly one Solid."
            continue
        }

        set solid_id [lindex $solid_ids 0]

        puts ""
        puts "------------------------------------------------------------"
        puts "Processing Solid ID = $solid_id"


        # ----------------------------------------------------------
        # Calculate geometric centroid
        # ----------------------------------------------------------

        if {[catch {
            set centroid [hm_getcentroid solids 1]
        } centroid_error]} {

            puts "ERROR: Cannot calculate centroid."
            puts $centroid_error

            incr failed_count

            if {![CYLAX::AskContinue $solid_id]} {
                break
            }

            continue
        }

        if {[llength $centroid] < 3} {

            puts "ERROR: Invalid centroid data."

            incr failed_count

            if {![CYLAX::AskContinue $solid_id]} {
                break
            }

            continue
        }


        # ----------------------------------------------------------
        # Define cylinder axis
        # ----------------------------------------------------------

        puts "Define the cylinder axis in the direction widget."
        puts "For N1-N2, select two points/nodes on the cylinder axis."
        puts "With projection_mode=closest, arrow direction is not important."

        if {[catch {
            set raw_direction [hm_getdirectionpanel \
                "Define cylinder axis" \
                N1N2N3 \
                0]
        } direction_error]} {

            puts "Direction definition was cancelled or failed."
            puts $direction_error

            incr failed_count

            if {![CYLAX::AskContinue $solid_id]} {
                break
            }

            continue
        }

        if {[catch {
            set direction \
                [CYLAX::ParseDirection $raw_direction]
        } parse_error]} {

            puts "ERROR: Invalid cylinder direction."
            puts $parse_error

            incr failed_count

            if {![CYLAX::AskContinue $solid_id]} {
                break
            }

            continue
        }


        # ----------------------------------------------------------
        # Calculate sufficient axis-line length
        # ----------------------------------------------------------

        set half_length [CYLAX::GetAxisHalfLength \
            $centroid \
            $bbox \
            $cfg(ray_length_scale) \
            $cfg(manual_ray_half_length)]

        puts "Axis line half-length = $half_length"


        # ----------------------------------------------------------
        # Find surface intersection
        # ----------------------------------------------------------

        if {[catch {
            set projection [CYLAX::FindProjection \
                $centroid \
                $direction \
                $half_length \
                $surface_ids \
                $cfg(projection_mode) \
                $cfg(minimum_projection_dist) \
                $cfg(allow_near_miss_fallback) \
                $cfg(near_miss_tolerance) \
                $cfg(keep_failed_axis_line)]
        } projection_error]} {

            puts "ERROR: Projection query failed."
            puts $projection_error

            incr failed_count

            if {![CYLAX::AskContinue $solid_id]} {
                break
            }

            continue
        }


        set status [dict get $projection status]


        # ----------------------------------------------------------
        # Failure diagnostics
        # ----------------------------------------------------------

        if {$status eq "NO_HIT"} {

            set debug_line_id \
                [dict get $projection line_id]

            set nearest_surface \
                [dict get $projection nearest_miss_surface]

            set nearest_distance \
                [dict get $projection nearest_miss_distance]

            puts ""
            puts "NO_HIT for Solid $solid_id"
            puts "Cylinder centroid = $centroid"
            puts "Cylinder axis     = $direction"

            if {
                $nearest_surface >= 0 &&
                $nearest_distance < 1.0e299
            } {
                puts "Nearest selected Surface ID = $nearest_surface"
                puts "Minimum line/surface distance = $nearest_distance"
            } else {
                puts "No usable surface-distance result was returned."
            }

            if {$debug_line_id ne ""} {
                puts "Failed axis line was retained."
                puts "Debug Line ID = $debug_line_id"
            }

            set centroid_local \
                [CYLAX::ToLocal \
                    $cfg(system_id) $centroid]

            set centroid_global_f \
                [CYLAX::FmtPoint \
                    $centroid $cfg(decimal_places)]

            set centroid_local_f \
                [CYLAX::FmtPoint \
                    $centroid_local $cfg(decimal_places)]

            set direction_f \
                [CYLAX::FmtPoint \
                    $direction $cfg(decimal_places)]

            lappend results [list \
                $solid_id \
                "NO_HIT" \
                {*}$centroid_global_f \
                {*}$centroid_local_f \
                {*}$direction_f \
                "" \
                $nearest_surface \
                "" "" "" \
                "" "" "" \
                "" \
                "" \
                $nearest_distance \
                "" \
                $debug_line_id]

            incr failed_count

            if {![CYLAX::AskContinue $solid_id]} {
                break
            }

            continue
        }


        # ----------------------------------------------------------
        # Successful or approximate result
        # ----------------------------------------------------------

        set surface_id \
            [dict get $projection surface_id]

        set projected_point \
            [dict get $projection point]

        set signed_distance \
            [dict get $projection signed_distance]

        set miss_distance \
            [dict get $projection miss_distance]

        if {$signed_distance >= 0.0} {
            set projection_side "positive"
        } else {
            set projection_side "negative"
        }

        set absolute_distance \
            [expr {abs($signed_distance)}]


        # ----------------------------------------------------------
        # Convert coordinates into System ID 1
        # ----------------------------------------------------------

        if {[catch {
            set centroid_local \
                [CYLAX::ToLocal \
                    $cfg(system_id) $centroid]

            set projected_local \
                [CYLAX::ToLocal \
                    $cfg(system_id) $projected_point]
        } local_error]} {

            puts "ERROR: Local coordinate conversion failed."
            puts $local_error

            incr failed_count

            if {![CYLAX::AskContinue $solid_id]} {
                break
            }

            continue
        }


        # ----------------------------------------------------------
        # Create projected node
        # ----------------------------------------------------------

        set created_node_id ""

        if {$cfg(create_projection_node)} {

            set px [lindex $projected_point 0]
            set py [lindex $projected_point 1]
            set pz [lindex $projected_point 2]

            if {[catch {
                *createnode $px $py $pz
                set created_node_id \
                    [hm_latestentityid nodes]
            } node_error]} {

                puts "WARNING: Projection point was found,"
                puts "but the output node could not be created."
                puts $node_error

                set created_node_id ""
            }
        }


        # ----------------------------------------------------------
        # Format output
        # ----------------------------------------------------------

        set centroid_global_f \
            [CYLAX::FmtPoint \
                $centroid $cfg(decimal_places)]

        set centroid_local_f \
            [CYLAX::FmtPoint \
                $centroid_local $cfg(decimal_places)]

        set direction_f \
            [CYLAX::FmtPoint \
                $direction $cfg(decimal_places)]

        set projected_global_f \
            [CYLAX::FmtPoint \
                $projected_point $cfg(decimal_places)]

        set projected_local_f \
            [CYLAX::FmtPoint \
                $projected_local $cfg(decimal_places)]

        set signed_distance_f \
            [CYLAX::Fmt \
                $signed_distance $cfg(decimal_places)]

        set absolute_distance_f \
            [CYLAX::Fmt \
                $absolute_distance $cfg(decimal_places)]

        set miss_distance_f \
            [CYLAX::Fmt \
                $miss_distance $cfg(decimal_places)]


        # ----------------------------------------------------------
        # Save result
        # ----------------------------------------------------------

        lappend results [list \
            $solid_id \
            $status \
            {*}$centroid_global_f \
            {*}$centroid_local_f \
            {*}$direction_f \
            $projection_side \
            $surface_id \
            {*}$projected_global_f \
            {*}$projected_local_f \
            $signed_distance_f \
            $absolute_distance_f \
            $miss_distance_f \
            $created_node_id \
            ""]


        # ----------------------------------------------------------
        # Console output
        # ----------------------------------------------------------

        puts ""
        puts "Projection completed."
        puts "Status               = $status"
        puts "Solid ID             = $solid_id"
        puts "Centroid global      = $centroid_global_f"
        puts "Centroid local       = $centroid_local_f"
        puts "Cylinder axis        = $direction_f"
        puts "Projection side      = $projection_side"
        puts "Surface ID           = $surface_id"
        puts "Projected global     = $projected_global_f"
        puts "Projected local      = $projected_local_f"
        puts "Axis distance        = $absolute_distance_f"
        puts "Line/surface miss    = $miss_distance_f"

        if {$created_node_id ne ""} {
            puts "Created node ID      = $created_node_id"
        }

        if {$status eq "EXACT"} {
            incr success_count
        } else {
            incr approximate_count
        }


        if {![CYLAX::AskContinue $solid_id]} {
            break
        }
    }


    # --------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------

    puts ""
    puts "============================================================"
    puts "Cylinder projection finished"
    puts "Exact projections       = $success_count"
    puts "Approximate projections = $approximate_count"
    puts "Failed projections      = $failed_count"
    puts "============================================================"


    CYLAX::ExportCSV

    puts ""
    puts "Script finished."
}


# ======================================================================
# Run
# ======================================================================

CYLAX::Main