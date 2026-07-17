# ============================================================
# HyperMesh Tcl Script
#
# 功能：
#   在局部坐标系 SID=1 下，按给定局部 XYZ 坐标创建节点
#
# 方法：
#   1. 在 SID=1 中创建临时几何点
#   2. 在临时点位置创建节点
#   3. 删除临时几何点
#   4. 验证节点在 SID=1 下的局部坐标
#
# 注意：
#   - 不交换 X、Y、Z
#   - 不手动进行局部到全局的坐标计算
#   - 重复运行会在相同位置创建重合节点
# ============================================================


# ------------------------------------------------------------
# 局部坐标系 ID
# ------------------------------------------------------------
set sid 1


# ------------------------------------------------------------
# 局部坐标验证误差
# ------------------------------------------------------------
set tolerance 1.0e-5


# ------------------------------------------------------------
# 节点坐标
#
# 每行格式：
#   {SID1局部X SID1局部Y SID1局部Z}
# ------------------------------------------------------------
set coords {
    {9.438146 4 2}
    {9.555818 4 7.6}
    {7.173491 4 12.3}
    {2.591164 4 14.9}
    {-2.59116 4 14.9}
    {-7.17349 4 12.3}
    {-9.55582 4 7.6}
    {-9.43815 4 2}

    {9.438146 8 2.5}
    {9.355818 8 7.9}
    {6.773491 8 11.8}
    {2.491164 8 13.9}
    {-2.49116 8 13.9}
    {-6.77349 8 11.8}
    {-9.35582 8 7.9}
    {-9.43815 8 2.5}

    {9.438146 12 3}
    {8.455818 12 8}
    {6.173491 12 10.7}
    {2.391164 12 12.5}
    {-2.39116 12 12.5}
    {-6.17349 12 10.7}
    {-8.45582 12 8}
    {-9.43815 12 3}

    {8.355818 16 4.7}
    {5.873491 16 8.3}
    {2.191164 16 10.2}
    {-2.19116 16 10.2}
    {-5.87349 16 8.3}
    {-8.35582 16 4.7}

    {6.855818 19 4}
    {4.773491 19 6.5}
    {1.991164 19 7.9}
    {-1.99116 19 7.9}
    {-4.77349 19 6.5}
    {-6.85582 19 4}

    {4.273491 21.2 3.8}
    {1.691164 21.5 5}
    {-1.69116 21.5 5}
    {-4.27349 21.2 3.8}
}


# ------------------------------------------------------------
# 检查局部坐标系 SID=1 是否存在
# ------------------------------------------------------------
if {[catch {
    set allSystemIDs [hm_entitylist systems id all]
} systemListError]} {
    error "Unable to read coordinate-system IDs: $systemListError"
}

if {[lsearch -exact $allSystemIDs $sid] < 0} {
    error "Coordinate system SID=$sid does not exist. Existing system IDs: $allSystemIDs"
}


# ------------------------------------------------------------
# 检查坐标系类型
#
# Type 0 = rectangular
# 当前数据为普通 XYZ 笛卡尔坐标，因此要求矩形坐标系
# ------------------------------------------------------------
set systemType ""

if {[catch {
    set systemType \
        [hm_getvalue systems id=$sid dataname=type]
} systemTypeError]} {
    error "Unable to read the type of SID=$sid: $systemTypeError"
}

if {$systemType != 0} {
    error "SID=$sid is not a rectangular coordinate system. System type=$systemType"
}


# ------------------------------------------------------------
# 输出坐标系基本信息
# ------------------------------------------------------------
set systemOrigin ""
set systemXAxis ""
set systemYAxis ""
set systemZAxis ""

catch {
    set systemOrigin \
        [hm_getvalue systems id=$sid dataname=origin]
}

catch {
    set systemXAxis \
        [hm_getvalue systems id=$sid dataname=xaxis]
}

catch {
    set systemYAxis \
        [hm_getvalue systems id=$sid dataname=yaxis]
}

catch {
    set systemZAxis \
        [hm_getvalue systems id=$sid dataname=zaxis]
}


puts ""
puts "================================================"
puts "Coordinate system confirmed"
puts "SID       = $sid"
puts "Type      = $systemType"
puts "Origin    = {$systemOrigin}"
puts "Local X   = {$systemXAxis}"
puts "Local Y   = {$systemYAxis}"
puts "Local Z   = {$systemZAxis}"
puts "Point Qty = [llength $coords]"
puts "================================================"


# ------------------------------------------------------------
# 初始化
# ------------------------------------------------------------
set createdNodeIDs {}
set failedRows {}
set row 0


# ------------------------------------------------------------
# 逐个创建节点
# ------------------------------------------------------------
foreach point $coords {

    incr row

    # --------------------------------------------------------
    # 检查当前坐标行
    # --------------------------------------------------------
    if {[llength $point] != 3} {
        lappend failedRows $row
        puts "ERROR: Row $row does not contain three coordinates: {$point}"
        continue
    }

    set localX [lindex $point 0]
    set localY [lindex $point 1]
    set localZ [lindex $point 2]

    if {![string is double -strict $localX] ||
        ![string is double -strict $localY] ||
        ![string is double -strict $localZ]} {

        lappend failedRows $row
        puts "ERROR: Invalid coordinate at row $row: {$point}"
        continue
    }


    # --------------------------------------------------------
    # 第一步：
    # 在 SID=1 的局部坐标位置创建临时几何点
    #
    # *createpoint x y z system_id
    # --------------------------------------------------------
    if {[catch {
        *createpoint $localX $localY $localZ $sid
    } pointCreateError]} {

        lappend failedRows $row

        puts "ERROR: Unable to create temporary point at row $row"
        puts "       Local XYZ = {$localX $localY $localZ}"
        puts "       Error     = $pointCreateError"

        continue
    }


    # 获取刚创建的临时点 ID
    set pointID [hm_latestentityid points]

    if {$pointID <= 0} {
        lappend failedRows $row
        puts "ERROR: Unable to obtain temporary point ID at row $row"
        continue
    }


    # --------------------------------------------------------
    # 第二步：
    # 在临时几何点的位置创建自由节点
    # --------------------------------------------------------
    *createmark points 1 $pointID

    if {[catch {
        *nodecreateatpointmark 1
    } nodeCreateError]} {

        lappend failedRows $row

        puts "ERROR: Unable to create node at row $row"
        puts "       Local XYZ = {$localX $localY $localZ}"
        puts "       Point ID  = $pointID"
        puts "       Error     = $nodeCreateError"

        # 删除失败行对应的临时点
        *createmark points 1 $pointID
        catch {
            *deletemark points 1
        }

        continue
    }


    # 获取刚创建的节点 ID
    set nodeID [hm_latestentityid nodes]

    if {$nodeID <= 0} {

        lappend failedRows $row

        puts "ERROR: Unable to obtain node ID at row $row"

        *createmark points 1 $pointID
        catch {
            *deletemark points 1
        }

        continue
    }


    # --------------------------------------------------------
    # 第三步：
    # 设置节点输入和输出坐标系为 SID=1
    #
    # 此操作不用于定位节点；
    # 节点位置已经由临时点确定
    # --------------------------------------------------------
    if {[catch {
        *setvalue nodes id=$nodeID \
            inputsystemid=$sid \
            outputsystemid=$sid
    } setSystemError]} {

        puts "WARNING: Node $nodeID was created at the correct position,"
        puts "         but assigning input/output SID failed."
        puts "         Error: $setSystemError"
    }


    # --------------------------------------------------------
    # 第四步：
    # 删除临时几何点，只留下节点
    # --------------------------------------------------------
    *createmark points 1 $pointID

    if {[catch {
        *deletemark points 1
    } pointDeleteError]} {

        puts "WARNING: Temporary point $pointID could not be deleted."
        puts "         Error: $pointDeleteError"
    }


    # --------------------------------------------------------
    # 第五步：
    # 将节点实际位置换算回 SID=1，验证是否正确
    # --------------------------------------------------------
    if {[catch {
        set checkedLocal \
            [hm_xformnodetolocal $nodeID $sid]
    } checkError]} {

        lappend failedRows $row

        puts "ERROR: Unable to verify node $nodeID"
        puts "       Error: $checkError"

        continue
    }


    set checkedX [lindex $checkedLocal 0]
    set checkedY [lindex $checkedLocal 1]
    set checkedZ [lindex $checkedLocal 2]


    set errorX [expr {abs($checkedX - $localX)}]
    set errorY [expr {abs($checkedY - $localY)}]
    set errorZ [expr {abs($checkedZ - $localZ)}]


    # --------------------------------------------------------
    # 判断误差
    # --------------------------------------------------------
    if {$errorX > $tolerance ||
        $errorY > $tolerance ||
        $errorZ > $tolerance} {

        lappend failedRows $row

        puts "ERROR: Coordinate verification failed at row $row"
        puts "       Node ID         = $nodeID"
        puts "       Requested local = {$localX $localY $localZ}"
        puts "       Checked local   = {$checkedLocal}"
        puts "       Error XYZ       = {$errorX $errorY $errorZ}"

        continue
    }


    # 保存成功创建的节点 ID
    lappend createdNodeIDs $nodeID


    puts "Node $nodeID created successfully:"
    puts "  Row             = $row"
    puts "  Requested local = {$localX $localY $localZ}"
    puts "  Checked local   = {$checkedLocal}"
    puts "  Error XYZ       = {$errorX $errorY $errorZ}"
}


# ------------------------------------------------------------
# 显示并定位刚创建的节点
# ------------------------------------------------------------
if {[llength $createdNodeIDs] > 0} {

    # 将本次创建的节点放入 Mark 1
    eval *createmark nodes 1 $createdNodeIDs

    # 打开实体显示
    catch {
        *showall
    }

    # 显示节点编号
    # 即使节点圆点很小，也可以通过编号看到节点位置
    catch {
        *numbersmark nodes 1 1
    }

    # 设置等轴测视图
    catch {
        *view iso1
    }

    # 将视图缩放到本次创建的节点
    catch {
        *window_entitymark nodes 1
    }

    # 用高亮方式检查节点
    catch {
        *reviewentitybymark 1 3 1 0
    }

    puts "Displayed node IDs: $createdNodeIDs"
}


# ------------------------------------------------------------
# 完成信息
# ------------------------------------------------------------
puts ""
puts "================================================"
puts "Node creation completed"
puts "Coordinate system SID = $sid"
puts "Requested nodes       = [llength $coords]"
puts "Created nodes         = [llength $createdNodeIDs]"
puts "Failed rows           = $failedRows"
puts "Created node IDs      = $createdNodeIDs"
puts "================================================"


# ------------------------------------------------------------
# 如果有失败坐标，最后显示错误
# ------------------------------------------------------------
if {[llength $failedRows] > 0} {

    error "Some coordinate rows failed verification: $failedRows"
}