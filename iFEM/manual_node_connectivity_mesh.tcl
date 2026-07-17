# =============================================================================
# manual_node_connectivity_mesh.tcl
#
# HyperMesh Tcl：完全手动指定节点连接关系
#
# 核心原则：
#   - 脚本不判断曲面轮廓；
#   - 脚本不自动猜测哪些节点应该连接；
#   - 脚本不创建、移动或投影节点；
#   - 每个 QUAD4/TRIA3 的节点连接关系由用户明确填写。
#
# 推荐工作流程：
#   第一次：
#       set RUN_MODE "report"
#       运行脚本，查看显示节点的 ID 和 SID=1 局部坐标。
#
#   第二次：
#       根据节点编号填写 QUADS 和 TRIAS；
#       set RUN_MODE "mesh"
#       再次运行脚本，创建单元。
#
# QUAD4 节点顺序：
#
#       n4 -------- n3
#       |            |
#       |            |
#       n1 -------- n2
#
#   必须沿单元边界连续排列，不能写成交叉顺序。
#
# TRIA3 节点顺序：
#
#              n3
#             /  \
#            /    \
#           n1----n2
#
# 节点顺序决定壳单元法向。法向相反时，可设置 FLIP_NORMAL=1，
# 或手动反转各单元节点顺序。
# =============================================================================


# =============================================================================
# 用户设置
# =============================================================================

# "report"：仅打印节点编号和坐标，不创建单元
# "mesh"  ：根据 QUADS / TRIAS 创建单元
set RUN_MODE "report"

# 用于打印节点局部坐标的坐标系
set SID 1

# report 模式及节点使用检查的范围：
#   "displayed" = 当前显示节点
#   "all"       = 模型全部节点
set NODE_SOURCE_MODE "displayed"

# 壳单元法向：
#   0 = 按 QUADS/TRIAS 中写入的节点顺序
#   1 = 全部反转
set FLIP_NORMAL 0

# 最大允许单元边长：
#   0.0 = 不检查
#   正数 = 任何单元边超过此长度时，在创建单元前停止
#
# 建议先设为一个合理值，可以防止误写远距离节点。
set MAX_ALLOWED_EDGE_LENGTH 0.0

# 最小允许单元面积：
# 小于该值视为退化单元。
set MIN_ELEMENT_AREA 1.0e-10

# 是否要求当前 NODE_SOURCE_MODE 范围内的每个节点都出现在连接表中：
#   1 = 要求全部使用
#   0 = 允许存在未使用节点
set REQUIRE_ALL_SOURCE_NODES_USED 1

# 是否显示节点编号
set SHOW_NODE_NUMBERS 1

# 是否显示新建单元编号
set SHOW_ELEMENT_NUMBERS 0

# 是否自动缩放到新建网格
set FIT_NEW_MESH 1


# =============================================================================
# 手动连接表——你只需要主要编辑这里
# =============================================================================

# 每一行是一个 QUAD4：
#     {n1 n2 n3 n4}
#
# 示例：
# set QUADS {
#     {101 102 202 201}
#     {102 103 203 202}
# }
#
# 注意：
#   n1 -> n2 -> n3 -> n4 必须沿四边形边界连续排列。
set QUADS {
    # {101 102 202 201}
    # {102 103 203 202}
}


# 每一行是一个 TRIA3：
#     {n1 n2 n3}
#
# 示例：
# set TRIAS {
#     {103 104 203}
# }
set TRIAS {
    # {103 104 203}
}


# 可选：明确列出“必须参与网格”的节点 ID。
#
# 若此列表非空，脚本只检查这里列出的节点是否被使用；
# 若此列表为空，并且 REQUIRE_ALL_SOURCE_NODES_USED=1，
# 脚本检查当前 displayed/all 范围内的全部节点。
#
# 示例：
# set REQUIRED_NODE_IDS {
#     101 102 103 104
#     201 202 203
# }
set REQUIRED_NODE_IDS {
}


# =============================================================================
# 基础工具
# =============================================================================

proc AppendUnique {inputList value} {
    if {[lsearch -exact $inputList $value] < 0} {
        lappend inputList $value
    }

    return $inputList
}


proc GetNodeXYZ {nodeID} {
    return [list \
        [hm_getvalue nodes id=$nodeID dataname=x] \
        [hm_getvalue nodes id=$nodeID dataname=y] \
        [hm_getvalue nodes id=$nodeID dataname=z]]
}


proc Distance3 {pointA pointB} {
    set dx [expr {[lindex $pointA 0] - [lindex $pointB 0]}]
    set dy [expr {[lindex $pointA 1] - [lindex $pointB 1]}]
    set dz [expr {[lindex $pointA 2] - [lindex $pointB 2]}]

    return [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
}


proc VectorSubtract {pointA pointB} {
    return [list \
        [expr {[lindex $pointA 0] - [lindex $pointB 0]}] \
        [expr {[lindex $pointA 1] - [lindex $pointB 1]}] \
        [expr {[lindex $pointA 2] - [lindex $pointB 2]}]]
}


proc CrossProduct {vectorA vectorB} {
    set ax [lindex $vectorA 0]
    set ay [lindex $vectorA 1]
    set az [lindex $vectorA 2]

    set bx [lindex $vectorB 0]
    set by [lindex $vectorB 1]
    set bz [lindex $vectorB 2]

    return [list \
        [expr {$ay*$bz - $az*$by}] \
        [expr {$az*$bx - $ax*$bz}] \
        [expr {$ax*$by - $ay*$bx}]]
}


proc VectorMagnitude {vector} {
    set x [lindex $vector 0]
    set y [lindex $vector 1]
    set z [lindex $vector 2]

    return [expr {sqrt($x*$x + $y*$y + $z*$z)}]
}


proc TriangleArea {point1 point2 point3} {
    set vector12 [VectorSubtract $point2 $point1]
    set vector13 [VectorSubtract $point3 $point1]

    set cross [CrossProduct $vector12 $vector13]

    return [expr {0.5 * [VectorMagnitude $cross]}]
}


proc ElementAreaFromNodes {elementType nodeIDs} {
    set point1 [GetNodeXYZ [lindex $nodeIDs 0]]
    set point2 [GetNodeXYZ [lindex $nodeIDs 1]]
    set point3 [GetNodeXYZ [lindex $nodeIDs 2]]

    if {[string equal $elementType "tria"]} {
        return [TriangleArea $point1 $point2 $point3]
    }

    set point4 [GetNodeXYZ [lindex $nodeIDs 3]]

    # QUAD4 使用两个三角形面积之和。
    return [expr {
        [TriangleArea $point1 $point2 $point3] +
        [TriangleArea $point1 $point3 $point4]
    }]
}


proc GetElementEdges {elementType nodeIDs} {
    if {[string equal $elementType "tria"]} {
        return [list \
            [list [lindex $nodeIDs 0] [lindex $nodeIDs 1]] \
            [list [lindex $nodeIDs 1] [lindex $nodeIDs 2]] \
            [list [lindex $nodeIDs 2] [lindex $nodeIDs 0]]]
    }

    return [list \
        [list [lindex $nodeIDs 0] [lindex $nodeIDs 1]] \
        [list [lindex $nodeIDs 1] [lindex $nodeIDs 2]] \
        [list [lindex $nodeIDs 2] [lindex $nodeIDs 3]] \
        [list [lindex $nodeIDs 3] [lindex $nodeIDs 0]]]
}


proc GetSourceNodeIDs {mode} {
    if {[string equal -nocase $mode "displayed"]} {
        *createmark nodes 1 "displayed"
        return [hm_getmark nodes 1]
    }

    if {[string equal -nocase $mode "all"]} {
        return [hm_entitylist nodes id all]
    }

    error "NODE_SOURCE_MODE must be displayed or all."
}


proc CreateQuad {nodeIDs flipNormal} {
    if {$flipNormal} {
        set nodeIDs [lreverse $nodeIDs]
    }

    eval *createlist node 1 $nodeIDs
    *createelement 104 1 1 0

    set elementID [hm_latestentityid elems]

    if {$elementID <= 0} {
        error "Unable to obtain newly created QUAD4 ID."
    }

    return $elementID
}


proc CreateTria {nodeIDs flipNormal} {
    if {$flipNormal} {
        set nodeIDs [lreverse $nodeIDs]
    }

    eval *createlist node 1 $nodeIDs
    *createelement 103 1 1 0

    set elementID [hm_latestentityid elems]

    if {$elementID <= 0} {
        error "Unable to obtain newly created TRIA3 ID."
    }

    return $elementID
}


# =============================================================================
# 1. 检查坐标系
# =============================================================================

set systemIDs [hm_entitylist systems id all]

if {[lsearch -exact $systemIDs $SID] < 0} {
    error "Coordinate system SID=$SID does not exist. Existing systems: $systemIDs"
}


# =============================================================================
# 2. Report 模式：打印节点编号和坐标
# =============================================================================

set sourceNodeIDs [GetSourceNodeIDs $NODE_SOURCE_MODE]

if {[llength $sourceNodeIDs] == 0} {
    error "No nodes were found in NODE_SOURCE_MODE=$NODE_SOURCE_MODE."
}

if {[string equal -nocase $RUN_MODE "report"]} {
    puts ""
    puts "================================================================"
    puts "NODE CONNECTIVITY REPORT"
    puts "Node source mode : $NODE_SOURCE_MODE"
    puts "Node count       : [llength $sourceNodeIDs]"
    puts "Local SID        : $SID"
    puts "================================================================"
    puts "Format:"
    puts "nodeID | localX localY localZ | globalX globalY globalZ"
    puts "================================================================"

    foreach nodeID [lsort -integer $sourceNodeIDs] {
        set localXYZ [hm_xformnodetolocal $nodeID $SID]
        set globalXYZ [GetNodeXYZ $nodeID]

        puts "$nodeID | $localXYZ | $globalXYZ"
    }

    if {$SHOW_NODE_NUMBERS} {
        eval *createmark nodes 1 $sourceNodeIDs

        catch {
            *numbersmark nodes 1 1
        }

        catch {
            *window_entitymark nodes 1
        }
    }

    puts "================================================================"
    puts "Report complete."
    puts ""
    puts "Next:"
    puts "  1. Copy node IDs into QUADS and TRIAS."
    puts "  2. Change RUN_MODE from report to mesh."
    puts "  3. Run this script again."
    puts "================================================================"

    return
}


if {![string equal -nocase $RUN_MODE "mesh"]} {
    error "RUN_MODE must be report or mesh."
}


# =============================================================================
# 3. 整理用户定义的单元
# =============================================================================

set elementSpecs {}
set connectivityNodeIDs {}
set validationErrors {}
set seenElementKeys [dict create]
set quadNumber 0
set triaNumber 0


foreach nodeIDs $QUADS {
    incr quadNumber

    if {[llength $nodeIDs] != 4} {
        lappend validationErrors \
            "QUAD $quadNumber must contain exactly 4 node IDs: {$nodeIDs}"
        continue
    }

    set uniqueIDs {}

    foreach nodeID $nodeIDs {
        set uniqueIDs [AppendUnique $uniqueIDs $nodeID]
    }

    if {[llength $uniqueIDs] != 4} {
        lappend validationErrors \
            "QUAD $quadNumber contains repeated node IDs: {$nodeIDs}"
        continue
    }

    set sortedKey [lsort -integer $nodeIDs]
    set key "quad:[join $sortedKey ,]"

    if {[dict exists $seenElementKeys $key]} {
        lappend validationErrors \
            "QUAD $quadNumber duplicates another quad node set: {$nodeIDs}"
        continue
    }

    dict set seenElementKeys $key 1

    lappend elementSpecs [list quad $nodeIDs "QUAD $quadNumber"]

    foreach nodeID $nodeIDs {
        set connectivityNodeIDs \
            [AppendUnique $connectivityNodeIDs $nodeID]
    }
}


foreach nodeIDs $TRIAS {
    incr triaNumber

    if {[llength $nodeIDs] != 3} {
        lappend validationErrors \
            "TRIA $triaNumber must contain exactly 3 node IDs: {$nodeIDs}"
        continue
    }

    set uniqueIDs {}

    foreach nodeID $nodeIDs {
        set uniqueIDs [AppendUnique $uniqueIDs $nodeID]
    }

    if {[llength $uniqueIDs] != 3} {
        lappend validationErrors \
            "TRIA $triaNumber contains repeated node IDs: {$nodeIDs}"
        continue
    }

    set sortedKey [lsort -integer $nodeIDs]
    set key "tria:[join $sortedKey ,]"

    if {[dict exists $seenElementKeys $key]} {
        lappend validationErrors \
            "TRIA $triaNumber duplicates another tria node set: {$nodeIDs}"
        continue
    }

    dict set seenElementKeys $key 1

    lappend elementSpecs [list tria $nodeIDs "TRIA $triaNumber"]

    foreach nodeID $nodeIDs {
        set connectivityNodeIDs \
            [AppendUnique $connectivityNodeIDs $nodeID]
    }
}


if {[llength $elementSpecs] == 0} {
    error "QUADS and TRIAS are empty. Fill the manual connectivity tables first."
}


# =============================================================================
# 4. 检查节点是否存在
# =============================================================================

set allModelNodeIDs [hm_entitylist nodes id all]

foreach nodeID $connectivityNodeIDs {
    if {[lsearch -exact $allModelNodeIDs $nodeID] < 0} {
        lappend validationErrors \
            "Connectivity references missing node ID $nodeID."
    }
}


# =============================================================================
# 5. 检查所有要求节点是否被使用
# =============================================================================

if {[llength $REQUIRED_NODE_IDS] > 0} {
    set requiredCheckIDs $REQUIRED_NODE_IDS
} elseif {$REQUIRE_ALL_SOURCE_NODES_USED} {
    set requiredCheckIDs $sourceNodeIDs
} else {
    set requiredCheckIDs {}
}

set unusedRequiredNodeIDs {}

foreach nodeID $requiredCheckIDs {
    if {[lsearch -exact $connectivityNodeIDs $nodeID] < 0} {
        lappend unusedRequiredNodeIDs $nodeID
    }
}

if {[llength $unusedRequiredNodeIDs] > 0} {
    lappend validationErrors \
        "Required nodes not used by QUADS/TRIAS: $unusedRequiredNodeIDs"
}


# =============================================================================
# 6. 单元几何检查
# =============================================================================

foreach spec $elementSpecs {
    set elementType [lindex $spec 0]
    set nodeIDs [lindex $spec 1]
    set label [lindex $spec 2]

    # 只有节点全部存在时才能进行坐标检查。
    set allExist 1

    foreach nodeID $nodeIDs {
        if {[lsearch -exact $allModelNodeIDs $nodeID] < 0} {
            set allExist 0
            break
        }
    }

    if {!$allExist} {
        continue
    }

    set area [ElementAreaFromNodes $elementType $nodeIDs]

    if {$area <= $MIN_ELEMENT_AREA} {
        lappend validationErrors \
            "$label is degenerate or nearly zero-area. Nodes={$nodeIDs}, area=$area"
    }

    if {$MAX_ALLOWED_EDGE_LENGTH > 0.0} {
        foreach edge [GetElementEdges $elementType $nodeIDs] {
            set nodeA [lindex $edge 0]
            set nodeB [lindex $edge 1]

            set pointA [GetNodeXYZ $nodeA]
            set pointB [GetNodeXYZ $nodeB]

            set edgeLength [Distance3 $pointA $pointB]

            if {$edgeLength > $MAX_ALLOWED_EDGE_LENGTH} {
                lappend validationErrors \
                    "$label edge {$nodeA $nodeB} length=$edgeLength exceeds MAX_ALLOWED_EDGE_LENGTH=$MAX_ALLOWED_EDGE_LENGTH"
            }
        }
    }
}


# =============================================================================
# 7. 验证失败时，不创建任何单元
# =============================================================================

if {[llength $validationErrors] > 0} {
    puts ""
    puts "================================================================"
    puts "MESH NOT CREATED — CONNECTIVITY VALIDATION FAILED"
    puts "================================================================"

    foreach message $validationErrors {
        puts "ERROR: $message"
    }

    puts "================================================================"

    error "Manual connectivity contains errors. No elements were created."
}


# =============================================================================
# 8. 创建单元
# =============================================================================

set createdElementIDs {}

foreach spec $elementSpecs {
    set elementType [lindex $spec 0]
    set nodeIDs [lindex $spec 1]
    set label [lindex $spec 2]

    if {[string equal $elementType "quad"]} {
        set elementID [CreateQuad $nodeIDs $FLIP_NORMAL]
    } else {
        set elementID [CreateTria $nodeIDs $FLIP_NORMAL]
    }

    lappend createdElementIDs $elementID

    puts "$label -> element $elementID, nodes={$nodeIDs}"
}


# =============================================================================
# 9. 创建后再次检查实际使用的节点
# =============================================================================

set actuallyUsedNodeIDs {}

foreach elementID $createdElementIDs {
    set elementNodeIDs \
        [hm_getvalue elems id=$elementID dataname=nodes]

    foreach nodeID $elementNodeIDs {
        set actuallyUsedNodeIDs \
            [AppendUnique $actuallyUsedNodeIDs $nodeID]
    }
}

set postCreationUnusedRequired {}

foreach nodeID $requiredCheckIDs {
    if {[lsearch -exact $actuallyUsedNodeIDs $nodeID] < 0} {
        lappend postCreationUnusedRequired $nodeID
    }
}


# =============================================================================
# 10. 显示结果
# =============================================================================

if {[llength $connectivityNodeIDs] > 0} {
    eval *createmark nodes 1 $connectivityNodeIDs

    catch {
        *makepreservednodes 1
    }

    if {$SHOW_NODE_NUMBERS} {
        catch {
            *numbersmark nodes 1 1
        }
    }
}

if {[llength $createdElementIDs] > 0} {
    eval *createmark elems 1 $createdElementIDs

    if {$SHOW_ELEMENT_NUMBERS} {
        catch {
            *numbersmark elems 1 1
        }
    }

    catch {
        *showall
    }

    if {$FIT_NEW_MESH} {
        catch {
            *window_entitymark elems 1
        }
    }
}


# =============================================================================
# 11. 总结
# =============================================================================

puts ""
puts "================================================================"
puts "MANUAL MESH CREATION COMPLETE"
puts "Defined QUAD4 count         : [llength $QUADS]"
puts "Defined TRIA3 count         : [llength $TRIAS]"
puts "Created element count       : [llength $createdElementIDs]"
puts "Connectivity node count     : [llength $connectivityNodeIDs]"
puts "Required node count         : [llength $requiredCheckIDs]"
puts "Unused required nodes       : $postCreationUnusedRequired"
puts "Created element IDs         : $createdElementIDs"
puts "================================================================"

if {[llength $postCreationUnusedRequired] == 0} {
    puts "CHECK OK: Every required node is used."
} else {
    puts "WARNING: Some required nodes were not used:"
    puts "         $postCreationUnusedRequired"
}

puts ""
puts "No automatic node-to-node connections were generated."
puts "Every element follows exactly the QUADS/TRIAS tables."
puts "================================================================"
