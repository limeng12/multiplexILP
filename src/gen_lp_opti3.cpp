#include <Rcpp.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <cmath>
#include <sstream>
#include <utility>
#include <regex>
#include <string> 

using namespace Rcpp;
using namespace std;


// 在文件开头，所有include之后，添加全局变量映射
// ============ 全局变量映射 ============
std::unordered_map<std::string, int> var_to_id;
std::unordered_map<int, std::string> id_to_var;
int next_var_id = 1;  // 从1开始编号

// 获取或创建变量ID
int get_or_create_var_id(const std::string& var_name) {
  auto it = var_to_id.find(var_name);
  if (it != var_to_id.end()) {
    return it->second;
  }
  int new_id = next_var_id++;
  var_to_id[var_name] = new_id;
  id_to_var[new_id] = var_name;
  return new_id;
}

std::pair<std::vector<int>, std::vector<int>> parse_to_pair_with_map(const std::vector<std::string>& inputs) {
  Rcout << "进入 parse_to_pair_with_map 函数" << endl;
  Rcout << "inputs.size() = " << inputs.size() << endl;
  
  std::vector<int> V1, V2;
  
  // 清除之前的映射
  var_to_id.clear();
  id_to_var.clear();
  next_var_id = 1;
  
  // 第一步：解析所有约束，建立冲突关系
  std::unordered_map<int, std::unordered_set<int>> conflicts;
  std::set<int> all_nodes_set;
  
  int matched_count = 0;
  int line_num = 0;
  
  for (const auto& input : inputs) {
    line_num++;
    std::string line = input;
    
    // 移除所有空格
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
    
    // 查找 "+" 和 "<=" 的位置
    size_t plus_pos = line.find('+');
    size_t le_pos = line.find("<=");
    
    if (plus_pos != std::string::npos && le_pos != std::string::npos && plus_pos < le_pos) {
      matched_count++;
      
      // 提取变量名
      std::string var1 = line.substr(0, plus_pos);
      std::string var2 = line.substr(plus_pos + 1, le_pos - plus_pos - 1);
      
      // 获取或创建变量ID
      int id1 = get_or_create_var_id(var1);
      int id2 = get_or_create_var_id(var2);

      // 记录冲突关系
      conflicts[id1].insert(id2);
      conflicts[id2].insert(id1);
      
      // 记录所有节点
      all_nodes_set.insert(id1);
      all_nodes_set.insert(id2);
      
      if (line_num <= 10) {
        Rcout << "第" << line_num << "行: '" << input << "' -> ";
        Rcout << "冲突: " << id1 << " -- " << id2 << " (原变量: " << var1 << "-" << var2 << ")" << endl;
      }
    }
  }
  
  // 第二步：生成补图的边（没有冲突的节点对）
  std::vector<int> all_nodes(all_nodes_set.begin(), all_nodes_set.end());
  std::sort(all_nodes.begin(), all_nodes.end());
  
  Rcout << "总节点数: " << all_nodes.size() << endl;
  
  int edge_count = 0;
  for (size_t i = 0; i < all_nodes.size(); i++) {
    int node1 = all_nodes[i];
    
    for (size_t j = i + 1; j < all_nodes.size(); j++) {
      int node2 = all_nodes[j];
      
      // 检查是否有冲突
      auto it1 = conflicts.find(node1);
      bool has_conflict = false;
      
      if (it1 != conflicts.end()) {
        if (it1->second.find(node2) != it1->second.end()) {
          has_conflict = true;
        }
      }
      
      // 如果没有冲突，则在补图中添加边
      if (!has_conflict) {
        V1.push_back(node1);
        V2.push_back(node2);
        edge_count++;
        
        if (edge_count <= 10) {
          Rcout << "补图边: " << node1 << " -- " << node2 
                << " (原变量: " << id_to_var[node1] << "-" << id_to_var[node2] << ")" << endl;
        }
      }
    }
  }
  
  Rcout << "解析完成，匹配到 " << matched_count << "/" << inputs.size() << " 个约束，";
  Rcout << "生成 " << V1.size() << " 条补图边，发现 " << var_to_id.size() << " 个不同变量" << endl;
  
  return std::make_pair(V1, V2);
}


struct DebugConfig {
  bool enabled = false;
  string output_dir = "./debug/";
  bool save_each_round = true;
  bool save_independent_sets = true;
};

void saveIndependentSets(const vector<unordered_set<int>>& independent_sets, 
                         const string& filename) {
  ofstream file(filename);
  if (!file.is_open()) {
    Rcerr << "无法打开文件: " << filename << endl;
    return;
  }
  
  file << "Independent Sets (共" << independent_sets.size() << "个):\n";
  file << "========================================\n";
  
  for (size_t i = 0; i < independent_sets.size(); i++) {
    // 将unordered_set转换为vector并排序以便显示
    vector<int> sorted_set(independent_sets[i].begin(), independent_sets[i].end());
    sort(sorted_set.begin(), sorted_set.end());
    
    file << "Set " << i << " (大小: " << sorted_set.size() << "):\n";
    file << "  节点: ";
    
    size_t count = 0;
    for (int node : sorted_set) {
      file << "X" << node;
      if (++count < sorted_set.size()) {
        file << ", ";
      }
      if (count % 10 == 0 && count < sorted_set.size()) {
        file << "\n       ";
      }
    }
    file << "\n\n";
  }
  
  file.close();
  Rcout << "独立集合已保存到: " << filename << endl;
}

// ============ 数据结构定义 ============
unordered_set<int> getAllNodes(const vector<int>& from, const vector<int>& to) {
  unordered_set<int> all_nodes;
  all_nodes.reserve(from.size() + to.size());
  for (int node : from) all_nodes.insert(node);
  for (int node : to) all_nodes.insert(node);
  return all_nodes;
}

unordered_map<int, unordered_set<int>> buildAdjacency(const vector<int>& from, const vector<int>& to) {
  unordered_map<int, unordered_set<int>> adj;
  for (size_t i = 0; i < from.size(); i++) {
    adj[from[i]].insert(to[i]);
    adj[to[i]].insert(from[i]);
  }
  return adj;
}

void ensureDebugDir(const string& dir_path) {
  // 创建目录（在R中应该已经存在）
}

List find_independent_sets_fixed(vector<int> from, vector<int> to, 
                                 bool debug = false, string debug_dir = "./debug/") {
  
  if (debug) {
    ensureDebugDir(debug_dir);
  }
  
  unordered_map<int, unordered_set<int>> adj;
  for (size_t i = 0; i < from.size(); i++) {
    adj[from[i]].insert(to[i]);
    adj[to[i]].insert(from[i]);
  }
  
  unordered_set<int> all_nodes;
  all_nodes.reserve(from.size() + to.size());
  for (int node : from) all_nodes.insert(node);
  for (int node : to) all_nodes.insert(node);
  
  Rcout << "总节点数: " << all_nodes.size() << endl;
  
  vector<unordered_set<int>> independent_sets;
  vector<int> node_to_set(1000000, -1);
  
  for (int node : all_nodes) {
    bool assigned = false;
    
    for (size_t set_id = 0; set_id < independent_sets.size(); set_id++) {
      bool can_add = true;
      
      // 检查节点是否与集合中任何节点相邻
      auto it = adj.find(node);
      if (it != adj.end()) {
        const unordered_set<int>& neighbors = it->second;
        for (int other_node : independent_sets[set_id]) {
          if (neighbors.find(other_node) != neighbors.end()) {
            can_add = false;
            break;
          }
        }
      }
      
      if (can_add) {
        independent_sets[set_id].insert(node);
        node_to_set[node] = set_id;
        assigned = true;
        break;
      }
    }
    
    if (!assigned) {
      unordered_set<int> new_set;
      new_set.insert(node);
      independent_sets.push_back(new_set);
      node_to_set[node] = independent_sets.size() - 1;
    }
  }
  
  if (debug) {
    saveIndependentSets(independent_sets, debug_dir + "independent_sets.txt");
  }
  
  List result;
  for (size_t i = 0; i < independent_sets.size(); i++) {
    // 转换为vector并排序以便返回
    vector<int> set_vec(independent_sets[i].begin(), independent_sets[i].end());
    sort(set_vec.begin(), set_vec.end());
    result.push_back(set_vec, "set_" + to_string(i));
  }
  
  Rcout <<  "总的独立集的个数："<< independent_sets.size() << "\n";
  
  
  return result;
}

struct Constraint {
  unordered_set<int> main_vars;        // 主变量集合
  int coefficient;                     // 系数
  unordered_set<int> other_vars;       // 其他变量集合
  int set_id;                          // 所属independent set
  bool active;                         // 是否活跃
  vector<int> source_ids;              // 来源约束ID（用于追踪）
  int round_created;                   // 创建轮次
  
  int size() const {
    return main_vars.size() + other_vars.size();
  }
  
  unordered_set<int> getAllVars() const {
    unordered_set<int> all_vars = main_vars;
    all_vars.insert(other_vars.begin(), other_vars.end());
    return all_vars;
  }
  
  string toString() const {
    string s = "";
    bool first = true;
    
    // 为了可读性，将变量排序后输出
    vector<int> sorted_main(main_vars.begin(), main_vars.end());
    sort(sorted_main.begin(), sorted_main.end());
    vector<int> sorted_other(other_vars.begin(), other_vars.end());
    sort(sorted_other.begin(), sorted_other.end());
    
    // 添加主变量（每个前面都有系数）
    for (int v : sorted_main) {
      if (!first) s += " + ";
      
      // 检查是否有变量名映射
      auto it = id_to_var.find(v);
      if (!id_to_var.empty() && it != id_to_var.end()) {
        // 有映射，使用原始变量名
        s += to_string(coefficient) + " " + it->second;
      } else {
        // 无映射，使用X加数字
        s += to_string(coefficient) + " X" + to_string(v);
      }
      first = false;
    }
    
    // 添加其他变量（没有系数）
    for (int v : sorted_other) {
      if (!first) s += " + ";
      
      // 检查是否有变量名映射
      auto it = id_to_var.find(v);
      if (!id_to_var.empty() && it != id_to_var.end()) {
        // 有映射，使用原始变量名
        s += it->second;
      } else {
        // 无映射，使用X加数字
        s += "X" + to_string(v);
      }
      first = false;
    }
    
    s += " <= " + to_string(coefficient);
    
    // 添加来源信息（调试用）
    if (!source_ids.empty()) {
      s += " [来自: ";
      for (size_t i = 0; i < source_ids.size(); i++) {
        s += to_string(source_ids[i]);
        if (i < source_ids.size() - 1) s += ",";
      }
      s += "]";
    }
    
    return s;
  }
  
  
  // 转换为R列表
  List toList() const {
    vector<int> main_vec(main_vars.begin(), main_vars.end());
    vector<int> other_vec(other_vars.begin(), other_vars.end());
    sort(main_vec.begin(), main_vec.end());
    sort(other_vec.begin(), other_vec.end());
    
    return List::create(
      Named("main_vars") = main_vec,
      Named("coefficient") = coefficient,
      Named("other_vars") = other_vec,
      Named("set_id") = set_id,
      Named("active") = active,
      Named("total_vars") = size(),
      Named("string") = toString(),
      Named("source_ids") = source_ids,
      Named("round_created") = round_created
    );
  }
  
  // 检查是否有重叠的主变量
  bool hasOverlappingMainVars(const Constraint& other) const {
    for (int v : main_vars) {
      if (other.main_vars.find(v) != other.main_vars.end()) {
        return true;
      }
    }
    return false;
  }
  
  // 获取所有共同变量（包括主变量和其他变量）
  unordered_set<int> getAllCommonVars(const Constraint& other) const {
    unordered_set<int> common_vars;
    const unordered_set<int>& all_vars1 = getAllVars();
    const unordered_set<int>& all_vars2 = other.getAllVars();
    
    // 遍历较小的集合
    const unordered_set<int>& smaller = (all_vars1.size() < all_vars2.size()) ? all_vars1 : all_vars2;
    const unordered_set<int>& larger = (all_vars1.size() < all_vars2.size()) ? all_vars2 : all_vars1;
    
    for (int v : smaller) {
      if (larger.find(v) != larger.end()) {
        common_vars.insert(v);
      }
    }
    
    return common_vars;
  }
  
  // 获取独特变量
  unordered_set<int> getUniqueVars(const Constraint& other) const {
    unordered_set<int> unique_vars = getAllVars();
    const unordered_set<int>& all_vars2 = other.getAllVars();
    
    // 移除共同的变量
    for (int v : all_vars2) {
      unique_vars.erase(v);
    }
    
    return unique_vars;
  }
};

struct ConstraintMerge {
  int constraint1_idx;
  int constraint2_idx;
  int reduction;
  unordered_set<int> common_main_vars;     // 共同的主变量（必须是独立且不重叠的）
  unordered_set<int> common_other_vars;    // 共同的其他变量
  unordered_set<int> unique1_vars;         // 约束1的独特变量
  unordered_set<int> unique2_vars;         // 约束2的独特变量
  int new_coefficient;
  
  bool operator<(const ConstraintMerge& other) const {
    if (reduction != other.reduction) {
      return reduction < other.reduction;
    }
    return (common_main_vars.size() + common_other_vars.size()) < 
      (other.common_main_vars.size() + other.common_other_vars.size());
  }
};

// 修改保存约束的函数，显示来源信息
void saveConstraintsWithSources(const vector<Constraint>& constraints, 
                                const string& filename,
                                const string& title = "约束",
                                int round = 0) {
  ofstream file(filename);
  if (!file.is_open()) {
    Rcerr << "无法打开文件: " << filename << endl;
    return;
  }
  
  file << title << " (共" << constraints.size() << "个，第" << round << "轮):\n";
  file << "========================================\n";
  
  int total_vars = 0;
  for (size_t i = 0; i < constraints.size(); i++) {
    const auto& cons = constraints[i];
    file << "约束 " << (i+1) << ": " << cons.toString() << "\n";
    
    file << "  主变量: ";
    vector<int> sorted_main(cons.main_vars.begin(), cons.main_vars.end());
    sort(sorted_main.begin(), sorted_main.end());
    for (int v : sorted_main) {
      file << "X" << v << " ";
    }
    file << "\n";
    
    file << "  其他变量: ";
    vector<int> sorted_other(cons.other_vars.begin(), cons.other_vars.end());
    sort(sorted_other.begin(), sorted_other.end());
    for (int v : sorted_other) {
      file << "X" << v << " ";
    }
    file << "\n";
    
    file << "  变量总数: " << cons.size() << "\n";
    file << "  所属集合: " << cons.set_id << "\n";
    file << "  是否活跃: " << (cons.active ? "是" : "否") << "\n";
    
    if (!cons.source_ids.empty()) {
      file << "  来源约束ID: ";
      for (size_t j = 0; j < cons.source_ids.size(); j++) {
        file << cons.source_ids[j];
        if (j < cons.source_ids.size() - 1) file << ", ";
      }
      file << "\n";
    }
    
    file << "  创建轮次: " << cons.round_created << "\n";
    file << "----------------------------------------\n";
    
    total_vars += cons.size();
  }
  
  file << "\n统计信息:\n";
  file << "  约束总数: " << constraints.size() << "\n";
  file << "  变量总数: " << total_vars << "\n";
  file << "  平均变量数: " << (constraints.empty() ? 0 : (double)total_vars / constraints.size()) << "\n";
  
  file.close();
  Rcout << "约束已保存到: " << filename << endl;
}

// ============ 生成约束 ============

vector<Constraint> generate_constraints(const vector<int>& from, 
                                        const vector<int>& to,
                                        const vector<int>& node_to_set) {
  
  unordered_set<int> all_nodes;
  all_nodes.reserve(from.size() + to.size());
  for (int node : from) all_nodes.insert(node);
  for (int node : to) all_nodes.insert(node);
  
  unordered_map<int, unordered_set<int>> adj;
  for (size_t i = 0; i < from.size(); i++) {
    adj[from[i]].insert(to[i]);
    adj[to[i]].insert(from[i]);
  }
  
  vector<int> nodes(all_nodes.begin(), all_nodes.end());
  sort(nodes.begin(), nodes.end());
  
  // 生成所有不兼容对
  vector<pair<int, int>> incompatible_pairs;
  for (size_t i = 0; i < nodes.size(); i++) {
    for (size_t j = i+1; j < nodes.size(); j++) {
      if (adj.find(nodes[i]) == adj.end() || 
          adj.at(nodes[i]).find(nodes[j]) == adj.at(nodes[i]).end()) {
        incompatible_pairs.push_back(make_pair(nodes[i], nodes[j]));
      }
    }
  }
  
  // 按主变量分组
  unordered_map<int, unordered_set<int>> constraints_by_main_var;
  for (const auto& pair : incompatible_pairs) {
    int main_var = min(pair.first, pair.second);
    int other_var = max(pair.first, pair.second);
    constraints_by_main_var[main_var].insert(other_var);
  }
  
  // 生成约束
  vector<Constraint> constraints;
  for (const auto& pair : constraints_by_main_var) {
    int main_var = pair.first;
    const unordered_set<int>& other_vars = pair.second;
    
    if (other_vars.empty()) continue;
    
    Constraint cons;
    cons.main_vars.insert(main_var);
    cons.coefficient = other_vars.size();
    cons.other_vars = other_vars;
    cons.active = true;
    
    if (main_var < (int)node_to_set.size()) {
      cons.set_id = node_to_set[main_var];
    } else {
      cons.set_id = -1;
    }
    
    constraints.push_back(cons);
  }
  
  return constraints;
}

// 计算约束列表的总变量数
int calculate_total_vars(const vector<Constraint>& constraints) {
  int total = 0;
  for (const auto& cons : constraints) {
    total += cons.size();
  }
  return total;
}

// ============ 改进的合并收益计算（修复丢失变量关系问题）- 使用unordered_set ============
ConstraintMerge calculate_merge_gain_corrected_fixed(const Constraint& c1, 
                                                     const Constraint& c2,
                                                     const vector<int>& node_to_set,
                                                     const unordered_map<int, unordered_set<int>>& adj) {
  ConstraintMerge merge;
  merge.constraint1_idx = -1;
  merge.constraint2_idx = -1;
  merge.reduction = -1;
  
  // 检查：两个约束的主变量不能相同
  for (int v : c1.main_vars) {
    if (c2.main_vars.find(v) != c2.main_vars.end()) {
      return merge;  // 有重叠的主变量，不能合并
    }
  }
  
  // 检查：所有主变量是否属于同一个independent set
  unordered_set<int> all_main_vars = c1.main_vars;
  all_main_vars.insert(c2.main_vars.begin(), c2.main_vars.end());
  
  int check_set_id = -1;
  for (int main_var : all_main_vars) {
    if (main_var >= (int)node_to_set.size()) continue;
    int current_set = node_to_set[main_var];
    
    if (check_set_id == -1) {
      check_set_id = current_set;
    } else if (current_set != check_set_id && current_set != -1) {
      return merge;
    }
  }
  
  // 获取c1和c2的所有other_vars
  const unordered_set<int>& other1 = c1.other_vars;
  const unordered_set<int>& other2 = c2.other_vars;
  
  // 计算共同的other_vars
  unordered_set<int> common_other;
  // 遍历较小的集合以获得更好的性能
  const unordered_set<int>& smaller = (other1.size() < other2.size()) ? other1 : other2;
  const unordered_set<int>& larger = (other1.size() < other2.size()) ? other2 : other1;
  
  for (int v : smaller) {
    if (larger.find(v) != larger.end()) {
      common_other.insert(v);
    }
  }
  
  // 如果共同other_vars太少，合并收益不大
  if (common_other.size() < 2) {
    return merge;
  }
  
  // 设置合并信息
  merge.common_main_vars = all_main_vars;
  merge.common_other_vars = common_other;
  
  // 新系数：共同other_vars的数量
  merge.new_coefficient = common_other.size();
  
  // 计算独特变量
  unordered_set<int> unique1_other = other1;
  unordered_set<int> unique2_other = other2;
  
  // 移除共同变量
  for (int v : common_other) {
    unique1_other.erase(v);
    unique2_other.erase(v);
  }
  
  // 从独特变量中移除已经在共同变量中的变量
  for (int v : all_main_vars) {
    unique1_other.erase(v);
    unique2_other.erase(v);
  }
  
  merge.unique1_vars = unique1_other;
  merge.unique2_vars = unique2_other;
  
  // 计算收益
  int original_vars = c1.size() + c2.size();
  
  // 合并后的情况：
  int merged_vars = merge.common_main_vars.size() + merge.common_other_vars.size();
  merged_vars += c1.main_vars.size() + merge.unique1_vars.size();
  merged_vars += c2.main_vars.size() + merge.unique2_vars.size();
  
  merge.reduction = original_vars - merged_vars;
  
  // 只考虑有显著收益的合并
  if (merge.reduction <= 2) {
    merge.reduction = -1;
  }
  
  return merge;
}

// 修改merge_two_constraints函数，保留来源信息
vector<Constraint> merge_two_constraints_fixed(const Constraint& c1,
                                               const Constraint& c2,
                                               const ConstraintMerge& merge_info,
                                               int round,
                                               const unordered_map<int, unordered_set<int>>& adj) {
  vector<Constraint> result;
  
  // 1. 创建合并后的新约束
  Constraint merged_constraint;
  merged_constraint.main_vars = merge_info.common_main_vars;
  merged_constraint.other_vars = merge_info.common_other_vars;
  merged_constraint.coefficient = merge_info.common_other_vars.size();
  merged_constraint.set_id = c1.set_id;
  merged_constraint.active = true;
  
  // 合并来源ID
  merged_constraint.source_ids = c1.source_ids;
  merged_constraint.source_ids.insert(merged_constraint.source_ids.end(), 
                                      c2.source_ids.begin(), c2.source_ids.end());
  merged_constraint.round_created = round;
  
  result.push_back(merged_constraint);
  
  // 2. 为c1的独特other_vars创建剩余约束
  if (!merge_info.unique1_vars.empty()) {
    Constraint remaining1;
    remaining1.main_vars = c1.main_vars;
    remaining1.other_vars = merge_info.unique1_vars;
    remaining1.coefficient = merge_info.unique1_vars.size();
    remaining1.set_id = c1.set_id;
    remaining1.active = true;
    remaining1.source_ids = c1.source_ids;
    remaining1.round_created = round;
    
    result.push_back(remaining1);
  }
  
  // 3. 为c2的独特other_vars创建剩余约束
  if (!merge_info.unique2_vars.empty()) {
    Constraint remaining2;
    remaining2.main_vars = c2.main_vars;
    remaining2.other_vars = merge_info.unique2_vars;
    remaining2.coefficient = merge_info.unique2_vars.size();
    remaining2.set_id = c2.set_id;
    remaining2.active = true;
    remaining2.source_ids = c2.source_ids;
    remaining2.round_created = round;
    
    result.push_back(remaining2);
  }
  
  return result;
}

// ============ 多轮合并算法 ============
// 修改多轮合并函数，传递邻接关系
List multi_round_merge_corrected_fixed(vector<Constraint> constraints,
                                       const vector<int>& node_to_set,
                                       const unordered_map<int, unordered_set<int>>& adj,
                                       bool debug = false,
                                       string debug_dir = "./debug/") {
  
  bool changed = true;
  int round = 1;
  vector<int> round_constraints;
  vector<int> round_total_vars;
  vector<int> round_reductions;
  
  // 为初始约束设置源ID
  for (size_t i = 0; i < constraints.size(); i++) {
    constraints[i].source_ids.push_back(i);
    constraints[i].round_created = 0;
  }
  
  // 记录初始状态
  int initial_total_vars = calculate_total_vars(constraints);
  round_constraints.push_back(constraints.size());
  round_total_vars.push_back(initial_total_vars);
  round_reductions.push_back(0);
  
  if (debug) {
    saveConstraintsWithSources(constraints, debug_dir + "round_0.txt", "初始约束", 0);
  }
  
  while (changed && round <= 5) {
    changed = false;
    int prev_total_vars = calculate_total_vars(constraints);
    Rcout << "第 " << round << " 轮合并，当前约束数: " << constraints.size() 
          << "，总变量数: " << prev_total_vars << endl;
    
    // 计算所有可能的合并
    priority_queue<ConstraintMerge> merge_queue;
    
    for (size_t i = 0; i < constraints.size(); i++) {
      if (!constraints[i].active) continue;
      
      for (size_t j = i+1; j < constraints.size(); j++) {
        if (!constraints[j].active) continue;
        
        if (constraints[i].set_id != constraints[j].set_id) {
          continue;
        }
        
        ConstraintMerge merge1 = calculate_merge_gain_corrected_fixed(
          constraints[i], constraints[j], node_to_set, adj);
        
        if (merge1.reduction > 0) {
          merge1.constraint1_idx = i;
          merge1.constraint2_idx = j;
          merge_queue.push(merge1);
        }
      }
    }
    
    if (merge_queue.empty()) {
      Rcout << "无更多合并机会" << endl;
      break;
    }
    
    // 标记已合并的约束
    vector<bool> merged(constraints.size(), false);
    vector<Constraint> new_constraints;
    
    // 贪心合并
    while (!merge_queue.empty()) {
      ConstraintMerge best_merge = merge_queue.top();
      merge_queue.pop();
      
      int idx1 = best_merge.constraint1_idx;
      int idx2 = best_merge.constraint2_idx;
      
      if (merged[idx1] || merged[idx2]) continue;
      
      // 执行合并
      vector<Constraint> merged_constraints = merge_two_constraints_fixed(
        constraints[idx1], constraints[idx2], best_merge, round, adj);
      
      // 添加到新约束列表
      new_constraints.insert(new_constraints.end(), 
                             merged_constraints.begin(), 
                             merged_constraints.end());
      
      // 标记为已合并
      merged[idx1] = true;
      merged[idx2] = true;
      changed = true;
    }
    
    // 添加未合并的约束
    for (size_t i = 0; i < constraints.size(); i++) {
      if (!merged[i]) {
        new_constraints.push_back(constraints[i]);
      }
    }
    
    // 计算本轮减少的变量数
    constraints = new_constraints;
    int current_total_vars = calculate_total_vars(constraints);
    int round_reduction = prev_total_vars - current_total_vars;
    
    // 记录本轮统计
    round_constraints.push_back(constraints.size());
    round_total_vars.push_back(current_total_vars);
    round_reductions.push_back(round_reduction);
    
    Rcout << "第 " << round << " 轮减少变量: " << round_reduction << " 个" << endl;
    
    if (debug) {
      saveConstraintsWithSources(constraints, 
                                 debug_dir + "round_" + to_string(round) + ".txt",
                                 "第" + to_string(round) + "轮合并后", round);
    }
    
    round++;
  }
  
  int final_total_vars = calculate_total_vars(constraints);
  Rcout << "合并完成，最终约束数: " << constraints.size() 
        << "，最终总变量数: " << final_total_vars << endl;
  
  // 将约束转换为R列表
  List constraints_list;
  for (size_t i = 0; i < constraints.size(); i++) {
    constraints_list.push_back(constraints[i].toList());
  }
  
  // 返回结果
  return List::create(
    Named("constraints") = constraints_list,
    Named("round_constraints") = round_constraints,
    Named("round_total_vars") = round_total_vars,
    Named("round_reductions") = round_reductions
  );
}

// 修改主优化函数，传递邻接关系
List optimize_constraints_corrected_fixed(vector<int> from, vector<int> to,
                                          bool debug = false,
                                          string debug_dir = "./debug/") {
  Rcout << "构建邻接关系" << endl;
  
  // 构建邻接关系
  unordered_map<int, unordered_set<int>> adj;
  for (size_t i = 0; i < from.size(); i++) {
    adj[from[i]].insert(to[i]);
    adj[to[i]].insert(from[i]);
  }
  Rcout << "开始构建独立集" << endl;
  
  // 1. 获取independent sets
  List sets_result = find_independent_sets_fixed(from, to, debug, debug_dir);
  vector<unordered_set<int>> independent_sets;
  
  for (int i = 0; i < sets_result.size(); i++) {
    vector<int> set_nodes = as<vector<int>>(sets_result[i]);
    unordered_set<int> set_set(set_nodes.begin(), set_nodes.end());
    independent_sets.push_back(set_set);
  }
  
  // 2. 构建节点到set的映射
  vector<int> node_to_set(10000, -1);
  for (size_t set_id = 0; set_id < independent_sets.size(); set_id++) {
    for (int node : independent_sets[set_id]) {
      node_to_set[node] = set_id;
    }
  }
  Rcout << "生成初始约束" << endl;
  
  // 3. 生成初始约束
  vector<Constraint> constraints = generate_constraints(from, to, node_to_set);
  
  if (debug) {
    Rcout << "生成 " << constraints.size() << " 个初始约束" << endl;
  }
  
  // 4. 多轮合并
  List merge_result = multi_round_merge_corrected_fixed(
    constraints, node_to_set, adj, debug, debug_dir);
  
  List constraints_list = as<List>(merge_result["constraints"]);
  vector<int> round_constraints = as<vector<int>>(merge_result["round_constraints"]);
  vector<int> round_total_vars = as<vector<int>>(merge_result["round_total_vars"]);
  vector<int> round_reductions = as<vector<int>>(merge_result["round_reductions"]);
  
  // 5. 计算统计信息
  int initial_total_vars = calculate_total_vars(constraints);
  
  // 从constraints_list中提取最终约束字符串
  vector<string> constraint_strings;
  int final_total_vars = 0;
  
  for (int i = 0; i < constraints_list.size(); i++) {
    List cons = as<List>(constraints_list[i]);
    constraint_strings.push_back(as<string>(cons["string"]));
    final_total_vars += as<int>(cons["total_vars"]);
  }
  
  double reduction_percentage = initial_total_vars > 0 ? 
  (1.0 - (double)final_total_vars / initial_total_vars) * 100.0 : 0.0;
  
  // 6. 返回结果
  return List::create(
    Named("independent_sets") = sets_result,
    Named("initial_constraint_count") = (int)constraints.size(),
    Named("final_constraint_count") = constraints_list.size(),
    Named("initial_total_vars") = initial_total_vars,
    Named("final_total_vars") = final_total_vars,
    Named("reduction_percentage") = reduction_percentage,
    Named("round_constraints") = round_constraints,
    Named("round_total_vars") = round_total_vars,
    Named("round_reductions") = round_reductions,
    Named("constraints") = constraint_strings,
    Named("constraints_details") = constraints_list
  );
}


void generate_final_lp2(IntegerVector from_vec, IntegerVector to_vec, 
                        string filename = "clique_optimized.lp") {
  
  
  std::vector<int> from = Rcpp::as<std::vector<int>>(from_vec);
  std::vector<int> to = Rcpp::as<std::vector<int>>(to_vec);
  
  var_to_id.clear();
  id_to_var.clear();
  // 获取优化结果
  List result = optimize_constraints_corrected_fixed(from, to, true, "./");
  vector<string> constraints = as<vector<string>>(result["constraints"]);
  int initial_total_vars = as<int>(result["initial_total_vars"]);
  int final_total_vars = as<int>(result["final_total_vars"]);
  double reduction_percentage = as<double>(result["reduction_percentage"]);
  vector<int> round_constraints = as<vector<int>>(result["round_constraints"]);
  vector<int> round_total_vars = as<vector<int>>(result["round_total_vars"]);
  vector<int> round_reductions = as<vector<int>>(result["round_reductions"]);
  
  // 获取所有节点
  unordered_set<int> all_nodes = getAllNodes(from, to);
  
  // 生成LP文件
  ofstream file(filename);
  
  file << "Maximize\n";
  file << "obj: ";
  
  vector<int> nodes(all_nodes.begin(), all_nodes.end());
  sort(nodes.begin(), nodes.end());
  
  for (size_t i = 0; i < nodes.size(); i++) {
    file << "X" << nodes[i];
    if (i < nodes.size() - 1) {
      file << " + ";
      if ((i + 1) % 10 == 0) file << "\n";
    }
  }
  file << "\n\n";
  
  file << "Subject To\n";
  
  // 写入优化后的约束
  int constraint_num = 1;
  for (string& cons_str : constraints) {
    size_t pos = cons_str.find(" [来自:");
    if (pos != string::npos) {
      cons_str = cons_str.substr(0, pos);
    }
    
    file << "c" << constraint_num << ": " << cons_str << "\n";
    constraint_num++;
  }
  
  // 添加至少选择1个节点的约束
  file << "\nmin_size: ";
  for (size_t i = 0; i < nodes.size(); i++) {
    file << "X" << nodes[i];
    if (i < nodes.size() - 1) {
      file << " + ";
      if ((i + 1) % 15 == 0) file << "\n";
    }
  }
  file << " >= 1\n";
  
  file << "\nBinary\n";
  for (int node : nodes) {
    file << "X" << node << "\n";
  }
  
  file << "\nEnd\n";
  file.close();
  
  Rcout << "\n===========================================\n";
  Rcout << "优化完成!\n";
  Rcout << "文件: " << filename << "\n";
  Rcout << "独立集合数: " << as<List>(result["independent_sets"]).size() << "\n";
  Rcout << "原始约束数: " << as<int>(result["initial_constraint_count"]) << "\n";
  Rcout << "最终约束数: " << as<int>(result["final_constraint_count"]) << "\n";
  Rcout << "原始总变量数: " << initial_total_vars << "\n";
  Rcout << "最终总变量数: " << final_total_vars << "\n";
  Rcout << "变量减少: " << reduction_percentage << "%\n";
  
  // 输出每轮统计
  Rcout << "\n每轮合并统计:\n";
  Rcout << "轮次\t约束数\t总变量数\t减少变量数\n";
  for (size_t i = 0; i < round_constraints.size(); i++) {
    Rcout << i << "\t" << round_constraints[i] << "\t" 
          << round_total_vars[i] << "\t" 
          << round_reductions[i] << "\n";
  }
  
  Rcout << "===========================================\n";
}


// [[Rcpp::export]]
CharacterVector convert_constraints(CharacterVector input) {
  
  vector<string> origin_cons=Rcpp::as<std::vector<std::string>>(input);
  //vector<int> from=parse_to_pair(origin_cons).first;
  //vector<int> to=parse_to_pair(origin_cons).second;
  
  
  Rcout <<  "原始总变量："<< origin_cons.size()*2 << "\n";
  
  
  auto parsed = parse_to_pair_with_map(origin_cons);
  vector<int> from = parsed.first;
  vector<int> to = parsed.second;
  Rcout << "parsing 字符串完成!\n";
  //for(int i=0;i<from.size();i++){
  //  Rcout<<from.at(i)<<" "<<to.at(i)<<endl;
    
  //}
  
  
  // 获取优化结果
  List result = optimize_constraints_corrected_fixed(from, to, false, "./");
  vector<string> constraints = as<vector<string>>(result["constraints"]);
  int initial_total_vars = as<int>(result["initial_total_vars"]);
  int final_total_vars = as<int>(result["final_total_vars"]);
  double reduction_percentage = as<double>(result["reduction_percentage"]);
  vector<int> round_constraints = as<vector<int>>(result["round_constraints"]);
  vector<int> round_total_vars = as<vector<int>>(result["round_total_vars"]);
  vector<int> round_reductions = as<vector<int>>(result["round_reductions"]);
  
  
  Rcout << "\n===========================================\n";
  Rcout << "优化完成!\n";
  Rcout << "独立集合数: " << as<List>(result["independent_sets"]).size() << "\n";
  Rcout << "原始约束数: " << as<int>(result["initial_constraint_count"]) << "\n";
  Rcout << "最终约束数: " << as<int>(result["final_constraint_count"]) << "\n";
  Rcout << "原始总变量数: " << initial_total_vars << "\n";
  Rcout << "最终总变量数: " << final_total_vars << "\n";
  Rcout << "变量减少: " << reduction_percentage << "%\n";
  
  // 输出每轮统计
  Rcout << "\n每轮合并统计:\n";
  Rcout << "轮次\t约束数\t总变量数\t减少变量数\n";
  for (size_t i = 0; i < round_constraints.size(); i++) {
    Rcout << i << "\t" << round_constraints[i] << "\t" 
          << round_total_vars[i] << "\t" 
          << round_reductions[i] << "\n";
  }
  
  Rcout << "===========================================\n";
  
  
  vector<string> final_constraints;
  
  for (string& cons_str : constraints) {
    size_t pos = cons_str.find(" [来自:");
    if (pos != string::npos) {
      cons_str = cons_str.substr(0, pos);
    }
    
    final_constraints.push_back(cons_str);
    
    //file << "c" << constraint_num << ": " << cons_str << "\n";
    //constraint_num++;
  }
  
  
  
  return wrap(final_constraints);
}




// 辅助函数：打印每轮统计
void print_round_stats(vector<int> from, vector<int> to) {
  List result = optimize_constraints_corrected_fixed(from, to, false, "./");
  
  vector<int> round_constraints = as<vector<int>>(result["round_constraints"]);
  vector<int> round_total_vars = as<vector<int>>(result["round_total_vars"]);
  vector<int> round_reductions = as<vector<int>>(result["round_reductions"]);
  
  Rcout << "\n优化过程每轮统计:\n";
  Rcout << "轮次\t约束数\t总变量数\t减少变量数\n";
  for (size_t i = 0; i < round_constraints.size(); i++) {
    Rcout << i << "\t" << round_constraints[i] << "\t" 
          << round_total_vars[i] << "\t" 
          << round_reductions[i] << "\n";
  }
  
  Rcout << "\n汇总统计:\n";
  Rcout << "原始约束数: " << as<int>(result["initial_constraint_count"]) << "\n";
  Rcout << "最终约束数: " << as<int>(result["final_constraint_count"]) << "\n";
  Rcout << "原始总变量数: " << as<int>(result["initial_total_vars"]) << "\n";
  Rcout << "最终总变量数: " << as<int>(result["final_total_vars"]) << "\n";
  Rcout << "变量减少: " << as<double>(result["reduction_percentage"]) << "%\n";
}