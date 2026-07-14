# ================================================================
# 逐一点击投影节点，读取其在 System ID = 1 中的坐标
# 适用版本：HyperMesh 2025
#
# 操作：
#   1. 运行脚本
#   2. 逐个点击投影节点
#   3. 每点击一个节点，控制台立即显示坐标
#   4. 按 Esc 或取消选择结束
#   5. 最后可导出 CSV
#
# 输出：
#   点击顺序
#   Node ID
#   全局 X/Y/Z
#   System 1 局部 X/Y/Z
# ================================================================


# ---------------- 用户设置 ----------------

set system_id 1

# 是否在结束后导出 CSV
# 1 = 导出
# 0 = 不导出
set export_csv 1

# 是否忽略重复点击的节点
# 1 = 忽略重复节点
# 0 = 允许重复记录
set skip_duplicates 1

# 输出小数位数
set decimal_places 6

# ------------------------------------------


# ================================================================
# 数字格式化函数
# ================================================================

proc FormatCoordinate {value digits} {
    return [format "%.${digits}f" $value]
}


# ================================================================
# 检查 System ID = 1 是否存在
# ================================================================

if {[catch {
    set checked_system_id \
        [hm_getvalue systems id=$system_id dataname=id]
} system_error]} {

    puts ""
    puts "错误：无法读取 System ID = $system_id"
    puts $system_error
    return
}

if {$checked_system_id eq ""} {
    puts ""
    puts "错误：System ID = $system_id 不存在。"
    return
}


# ================================================================
# 初始化结果
# ================================================================

set result_rows {}

lappend result_rows [list \
    "Sequence" \
    "Node_ID" \
    "Global_X" \
    "Global_Y" \
    "Global_Z" \
    "System_1_X" \
    "System_1_Y" \
    "System_1_Z"]

set selected_node_ids {}
set sequence_number 0


puts ""
puts "=============================================================="
puts "逐一读取投影节点坐标"
puts "参考坐标系：System ID = $system_id"
puts "=============================================================="
puts ""
puts "请逐个点击投影节点。"
puts "每次只能点击一个节点。"
puts "按 Esc 或取消选择即可结束。"
puts ""


# ================================================================
# 循环逐一选择节点
# ================================================================

while {1} {

    # 单实体选择器
    *createentitypanel nodes \
        "点击一个投影节点；按 Esc 或取消可结束"

    # 必须紧接着读取刚选中的实体
    set node_id [hm_info lastselectedentity nodes]

    # 返回 0 表示取消或没有选择
    if {$node_id == 0} {
        puts ""
        puts "节点选择结束。"
        break
    }


    # ------------------------------------------------------------
    # 检查重复选择
    # ------------------------------------------------------------

    if {$skip_duplicates} {

        if {[lsearch -exact $selected_node_ids $node_id] >= 0} {
            puts ""
            puts "节点 $node_id 已经读取过，本次忽略。"
            continue
        }
    }


    # ------------------------------------------------------------
    # 获取全局坐标
    #
    # System ID = 0 代表全局坐标系
    # ------------------------------------------------------------

    if {[catch {
        set global_coordinates \
            [hm_xformnodetolocal $node_id 0]
    } global_error]} {

        puts ""
        puts "错误：无法读取 Node $node_id 的全局坐标。"
        puts $global_error
        continue
    }

    if {[llength $global_coordinates] < 3} {
        puts ""
        puts "错误：Node $node_id 的全局坐标数据不完整。"
        continue
    }


    # ------------------------------------------------------------
    # 获取 System ID = 1 中的局部坐标
    # ------------------------------------------------------------

    if {[catch {
        set local_coordinates \
            [hm_xformnodetolocal $node_id $system_id]
    } local_error]} {

        puts ""
        puts "错误：无法将 Node $node_id 转换到 System $system_id。"
        puts $local_error
        continue
    }

    if {[llength $local_coordinates] < 3} {
        puts ""
        puts "错误：Node $node_id 的局部坐标数据不完整。"
        continue
    }


    # ------------------------------------------------------------
    # 拆分坐标
    # ------------------------------------------------------------

    set global_x [lindex $global_coordinates 0]
    set global_y [lindex $global_coordinates 1]
    set global_z [lindex $global_coordinates 2]

    set local_x [lindex $local_coordinates 0]
    set local_y [lindex $local_coordinates 1]
    set local_z [lindex $local_coordinates 2]


    # ------------------------------------------------------------
    # 格式化
    # ------------------------------------------------------------

    set global_x_f \
        [FormatCoordinate $global_x $decimal_places]

    set global_y_f \
        [FormatCoordinate $global_y $decimal_places]

    set global_z_f \
        [FormatCoordinate $global_z $decimal_places]

    set local_x_f \
        [FormatCoordinate $local_x $decimal_places]

    set local_y_f \
        [FormatCoordinate $local_y $decimal_places]

    set local_z_f \
        [FormatCoordinate $local_z $decimal_places]


    # ------------------------------------------------------------
    # 保存结果
    # ------------------------------------------------------------

    incr sequence_number

    lappend selected_node_ids $node_id

    lappend result_rows [list \
        $sequence_number \
        $node_id \
        $global_x_f \
        $global_y_f \
        $global_z_f \
        $local_x_f \
        $local_y_f \
        $local_z_f]


    # ------------------------------------------------------------
    # 控制台立即输出
    # ------------------------------------------------------------

    puts ""
    puts "--------------------------------------------------------------"
    puts "点击顺序       = $sequence_number"
    puts "Node ID        = $node_id"
    puts "全局坐标       = $global_x_f  $global_y_f  $global_z_f"
    puts "System 1 坐标  = $local_x_f  $local_y_f  $local_z_f"
    puts "--------------------------------------------------------------"
}


# ================================================================
# 控制台汇总表
# ================================================================

puts ""
puts ""
puts "=============================================================================================="

puts [format \
    "%-10s %-12s %-15s %-15s %-15s %-15s %-15s %-15s" \
    "Sequence" \
    "Node ID" \
    "Global X" \
    "Global Y" \
    "Global Z" \
    "System 1 X" \
    "System 1 Y" \
    "System 1 Z"]

puts "----------------------------------------------------------------------------------------------"

foreach row [lrange $result_rows 1 end] {

    puts [format \
        "%-10s %-12s %-15s %-15s %-15s %-15s %-15s %-15s" \
        [lindex $row 0] \
        [lindex $row 1] \
        [lindex $row 2] \
        [lindex $row 3] \
        [lindex $row 4] \
        [lindex $row 5] \
        [lindex $row 6] \
        [lindex $row 7]]
}

puts "=============================================================================================="
puts "读取完成，共记录 $sequence_number 个节点。"
puts ""


# ================================================================
# 导出 CSV
# ================================================================

if {$export_csv && $sequence_number > 0} {

    set csv_path ""

    if {[catch {

        set csv_path [tk_getSaveFile \
            -title "保存投影节点坐标" \
            -defaultextension ".csv" \
            -initialfile "Projected_Nodes_System_1.csv" \
            -filetypes {
                {"CSV 文件" {.csv}}
                {"所有文件" {*}}
            }]

    } dialog_error]} {

        puts "无法打开文件保存窗口："
        puts $dialog_error
        set csv_path ""
    }


    if {$csv_path ne ""} {

        set file_handle ""

        if {[catch {

            set file_handle [open $csv_path "w"]

            fconfigure $file_handle \
                -encoding utf-8 \
                -translation lf

            # UTF-8 BOM，方便 Excel 正确显示中文
            puts -nonewline $file_handle "\ufeff"

            foreach row $result_rows {
                puts $file_handle [join $row ","]
            }

            close $file_handle
            set file_handle ""

        } write_error]} {

            if {$file_handle ne ""} {
                catch {close $file_handle}
            }

            puts ""
            puts "CSV 写入失败："
            puts $write_error

        } else {

            puts ""
            puts "CSV 已保存到："
            puts $csv_path
        }

    } else {

        puts "已取消 CSV 保存。"
    }
}

puts ""
puts "脚本执行结束。"