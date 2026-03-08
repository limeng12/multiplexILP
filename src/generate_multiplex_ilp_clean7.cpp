#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>
#include <cmath>
#include <queue>
#include <bitset>
#include <string> 

using namespace Rcpp;
using namespace std;


constexpr size_t MAX_PRIMERS = 3001;  // 不超过5000，设为5000足够

//[[Rcpp::export(rng = false)]]
List generate_multiplex_ilp_single_primers_direct(
    StringVector primers,               // 引物序列
    StringVector primer_ids,            // 引物ID
    StringVector primer53,              // "5"或"3"表示引物端
    NumericVector primer_temp,
    StringVector locus_ids,             // 引物所属位点
    NumericVector distances,            // 引物到locus的距离
    NumericVector copy_length_max,      // 位点最大可能长度
    NumericVector copy_length_min,      // 位点最小可能长度
    NumericVector repeat_distances,     // 重复区域的距离（每个位点一个值）
    StringVector compatible_pairs,      // 兼容的引物对（单个引物之间）
    NumericVector num_channels_v,                   // 通道数量
    NumericVector size_gap_v = 10                   // 大小间隔
) {
  
  int num_channels=num_channels_v.at(0);
  int size_gap =size_gap_v.at(0);
  
  int n = primers.size();  // 引物总数
  
  // 1. 建立引物信息结构
  struct PrimerInfo {
    string primer_id;
    string primer_seq;
    string locus;
    string end_type;  // "5" 或 "3"
    double distance;
    double copy_max;
    double copy_min;
    double temp;
    double repeat_dist;
    int index;
  };
  
  vector<PrimerInfo> all_primers;
  unordered_map<string, int> primer_id_to_idx;
  unordered_map<string, vector<int>> locus_to_primers;      // 位点到所有引物的映射
  unordered_map<string, vector<int>> locus_to_primers_5;    // 位点到5'端引物的映射
  unordered_map<string, vector<int>> locus_to_primers_3;    // 位点到3'端引物的映射
  unordered_map<string, int> locus_to_index;
  vector<string> loci;
  
  // 收集所有位点信息
  for(int i = 0; i < n; i++) {
    string locus = as<string>(locus_ids[i]);
    if(locus_to_index.find(locus) == locus_to_index.end()) {
      locus_to_index[locus] = loci.size();
      loci.push_back(locus);
    }
  }
  
  int num_loci = loci.size();
  
  // 处理每个引物
  for(int i = 0; i < n; i++) {
    PrimerInfo primer;
    primer.primer_id = as<string>(primer_ids[i]);
    primer.primer_seq = as<string>(primers[i]);
    primer.locus = as<string>(locus_ids[i]);
    primer.end_type = as<string>(primer53[i]);
    primer.distance = distances[i];
    primer.copy_max = copy_length_max[i];
    primer.copy_min = copy_length_min[i];
    primer.repeat_dist = repeat_distances[i];
    primer.temp=primer_temp[i];
    primer.index = all_primers.size();
    
    all_primers.push_back(primer);
    primer_id_to_idx[primer.primer_id] = primer.index;
    
    // 更新各种映射
    locus_to_primers[primer.locus].push_back(primer.index);
    
    if(primer.end_type == "5") {
      locus_to_primers_5[primer.locus].push_back(primer.index);
    } else if(primer.end_type == "3") {
      locus_to_primers_3[primer.locus].push_back(primer.index);
    }
  }
  
  int num_primers = all_primers.size();
  
  // 2. 建立引物之间的兼容性集合
  set<pair<int, int>> compatible_set;
  
  for(int i = 0; i < compatible_pairs.size(); i++) {
    string comp_str = as<string>(compatible_pairs[i]);
    size_t pos = comp_str.find("&");
    if(pos != string::npos) {
      string primer1 = comp_str.substr(0, pos);
      string primer2 = comp_str.substr(pos + 1);
      
      if(primer_id_to_idx.find(primer1) != primer_id_to_idx.end() &&
         primer_id_to_idx.find(primer2) != primer_id_to_idx.end()) {
        int idx1 = primer_id_to_idx[primer1];
        int idx2 = primer_id_to_idx[primer2];
        
        if( abs(all_primers[idx1].temp-all_primers[idx1].temp)>2.5) continue;
        
        //if(idx1 < idx2) {
        compatible_set.insert(make_pair(idx1, idx2));
        //} else {
        compatible_set.insert(make_pair(idx2, idx1));
        //}
      }
    }
  }
  
  // 3. 生成ILP模型
  string lp_filename = "single_primer_direct_model.lp";
  ofstream lp(lp_filename);
  lp << fixed << setprecision(2);
  
  // 3.1 目标函数：最大化总权重
  // 每个选中引物的系数 = 50000 - 距离×10
  lp << "Maximize\n";
  lp << "obj: ";
  
  // 变量定义
  vector<string> X_vars;       // X_i: 引物选择变量（二进制）
  vector<string> P_vars;       // P_k: 位点被选中的变量
  vector<string> C_vars;       // C_{k,t}: 位点k分配到通道t（二进制）
  vector<string> D_vars;       // D_{k,l}: 位点k和l在同一通道
  vector<string> E_vars;       // E_{k,l,t}: 辅助变量
  unordered_set<string> Z_vars; // Z_k: 位点k未被选中
  
  // 创建X_i变量
  for(int i = 0; i < num_primers; i++) {
    string X_var = "X" + to_string(i);
    X_vars.push_back(X_var);
  }
  
  // 目标函数：Σ[(50000 - distance_i×10) × X_i]
  bool first_term = true;
  for(int i = 0; i < num_primers; i++) {
    int coefficient = 50000 - all_primers[i].distance * 1.0;
    
    if(!first_term) {
      lp << " + ";
    }
    
    lp << coefficient << " " << "X" << i;
    first_term = false;
    
    // 每5个引物换行
    if((i + 1) % 5 == 0 && i != num_primers - 1) {
      lp << "\n      ";
    }
  }
  lp << "\n\n";
  
  // 约束部分
  lp << "Subject To\n";
  int constr_count = 0;
  
  // 约束1: 每个位点最多只能选择一个5'端引物
  for(const auto& locus_entry : locus_to_primers_5) {
    const vector<int>& primers_5 = locus_entry.second;
    
    if(primers_5.size() > 1) {
      for(size_t idx = 0; idx < primers_5.size(); idx++) {
        lp << "X" << primers_5[idx];
        if(idx != primers_5.size() - 1) lp << " + ";
      }
      lp << " <= 1\n";
      constr_count++;
    }
  }
  
  // 约束2: 每个位点最多只能选择一个3'端引物
  for(const auto& locus_entry : locus_to_primers_3) {
    const vector<int>& primers_3 = locus_entry.second;
    
    if(primers_3.size() > 1) {
      for(size_t idx = 0; idx < primers_3.size(); idx++) {
        lp << "X" << primers_3[idx];
        if(idx != primers_3.size() - 1) lp << " + ";
      }
      lp << " <= 1\n";
      constr_count++;
    }
  }
  
  // 创建P_k变量（位点是否被选中）
  vector<string> P_var_list(num_loci);
  for(int k = 0; k < num_loci; k++) {
    string P_var = "P" + to_string(k);
    P_var_list[k] = P_var;
    P_vars.push_back(P_var);
  }
  
  // 约束3: 定义P_k变量（位点被选中当且仅当同时选中了5'和3'引物）
  for(int k = 0; k < num_loci; k++) {
    string locus = loci[k];
    
    // 获取该位点的5'和3'引物
    if(locus_to_primers_5.find(locus) == locus_to_primers_5.end() ||
       locus_to_primers_3.find(locus) == locus_to_primers_3.end()) {
      continue;
    }
    
    const vector<int>& primers_5 = locus_to_primers_5[locus];
    const vector<int>& primers_3 = locus_to_primers_3[locus];
    
    if(primers_5.empty() || primers_3.empty()) {
      continue;
    }
    
    string P_var = P_var_list[k];
    
    
    //左侧引物数量等于右侧引物数量
    lp << "";
    bool first = true;
    
    // 第一部分：所有5'引物的和
    for(int idx_5 : primers_5) {
      if(!first) lp << " + ";
      lp << "X" << idx_5;
      first = false;
      constr_count++;
    }
    
    // 第二部分：减去所有3'引物
    for(int idx_3 : primers_3) {
      // 这里直接写负号，因为前面已经有项了
      lp << " - X" << idx_3;
      constr_count++;
    }
    
    lp << " == 0\n";
    
    
    
    
    
    // P_k ≤ X_i5 (对于每个5'引物)
    lp << P_var;
    for(int idx_5 : primers_5) {
      lp << " - X" << idx_5;
      constr_count++;
    }
    lp << " <= 0\n";
    
    
    // P_k ≤ X_i3 (对于每个3'引物)
    lp << P_var;
    for(int idx_3 : primers_3) {
      lp << " - X" << idx_3 ;
      constr_count++;
    }
    lp << " <= 0\n";
    
    
    // P_k ≥ X_i5 + X_i3 - 1 (对于所有5'+3'组合)
    // 这确保了只有当同时选中5'和3'引物时，P_k=1
    lp << P_var;
    for(int idx_5 : primers_5) {
      lp << " - X" << idx_5;
      constr_count++;
    }
    for(int idx_3 : primers_3) {
      lp << " - X" << idx_3;
      constr_count++;
    }
    
    lp << " >= -1\n";
    
    
  }
  
  // 创建C_{k,t}变量（位点-通道分配变量）
  vector<vector<string>> C_vars_matrix(num_loci, vector<string>(num_channels));
  for(int k = 0; k < num_loci; k++) {
    for(int t = 0; t < num_channels; t++) {
      string C_var = "C" + to_string(k) + "_" + to_string(t);
      C_vars_matrix[k][t] = C_var;
      C_vars.push_back(C_var);
    }
  }
  
  // 约束4: 每个位点最多只能分配到一个通道
  for(int k = 0; k < num_loci; k++) {
    for(int t = 0; t < num_channels; t++) {
      lp << C_vars_matrix[k][t];
      if(t != num_channels - 1) lp << " + ";
    }
    lp << " <= 1\n";
    constr_count++;
  }
  
  // 约束5: 如果位点k被选中，那么它必须被分配到一个通道
  for(int k = 0; k < num_loci; k++) {
    string P_var = P_var_list[k];
    
    lp << P_var;
    for(int t = 0; t < num_channels; t++) {
      lp << " - " << C_vars_matrix[k][t];
    }
    lp << " <= 0\n";
    constr_count++;
  }
  
  // 约束6: 如果位点k被分配到通道t，那么它必须被选中
  for(int k = 0; k < num_loci; k++) {
    for(int t = 0; t < num_channels; t++) {
      string C_var = C_vars_matrix[k][t];
      string P_var = P_var_list[k];
      
      lp << C_var << " - " << P_var << " <= 0\n";
      constr_count++;
    }
  }
  
  // 约束7: 定义位点的Z_k变量 (Z_k = 1 - P_k)
  for(int k = 0; k < num_loci; k++) {
    string Z_var = "Z" + to_string(k);
    Z_vars.insert(Z_var);
    string P_var = P_var_list[k];
    
    lp << Z_var << " + " << P_var << " = 1\n";
    constr_count++;
  }
  
  // 约束8: 定义位点之间的D_{k,l}变量
  // D_{k,l} = 1 当且仅当 存在某个通道t使得 C_{k,t} = C_{l,t} = 1
  for(int k = 0; k < num_loci; k++) {
    for(int l = k + 1; l < num_loci; l++) {
      string D_var = "D" + to_string(k) + "_" + to_string(l);
      D_vars.push_back(D_var);
      
      // 对于每个通道t，定义辅助变量E_{k,l,t} = C_{k,t} * C_{l,t}
      for(int t = 0; t < num_channels; t++) {
        string E_var = "E" + to_string(k) + "_" + to_string(l) + "_" + to_string(t);
        E_vars.push_back(E_var);
        string C_kt = C_vars_matrix[k][t];
        string C_lt = C_vars_matrix[l][t];
        
        // E_{k,l,t} ≤ C_{k,t}
        lp << E_var << " - " << C_kt << " <= 0\n";
        // E_{k,l,t} ≤ C_{l,t}
        lp << E_var << " - " << C_lt << " <= 0\n";
        // E_{k,l,t} ≥ C_{k,t} + C_{l,t} - 1
        lp << E_var << " - " << C_kt << " - " << C_lt << " >= -1\n";
        constr_count += 3;
        
        // D_{k,l} ≥ E_{k,l,t} 对于所有t
        lp << D_var << " - " << E_var << " >= 0\n";
        constr_count++;
      }
      
      // D_{k,l} ≤ ∑_t E_{k,l,t}
      lp << D_var;
      for(int t = 0; t < num_channels; t++) {
        string E_var = "E" + to_string(k) + "_" + to_string(l) + "_" + to_string(t);
        lp << " - " << E_var;
      }
      lp << " <= 0\n";
      constr_count++;
    }
  }
  
  // 约束9: 不兼容引物不能共存
  for(int i = 0; i < num_primers; i++) {
    
    string X_i = "X" + to_string(i);
    
    string concat_str="";
    
    int n_uncompitable=0;
    
    
    for(int j = i + 1; j < num_primers; j++) {
      // 同一引物位点的不同端引物是兼容的
      if(all_primers[i].locus == all_primers[j].locus) {
        // 但同一端的不允许同时选中
        if(all_primers[i].end_type == all_primers[j].end_type) {
          //lp << "X" << i << " + X" << j << " <= 1\n";
          
          concat_str+=" + X"+to_string(j);
          n_uncompitable++;
          constr_count++;
        }
        continue;
      }
      
      bool is_compatible = (compatible_set.find(make_pair(i, j)) != compatible_set.end());
      if(is_compatible) continue;
      
      //lp << "X" << i << " + X" << j << " <= 1\n";
      concat_str+=" + X"+to_string(j);
      n_uncompitable++;
      
      constr_count++;
    }
    
    if(n_uncompitable>0){
      lp << n_uncompitable <<" "<<X_i<< concat_str <<" <= "<<n_uncompitable <<"\n";
    }
    
  }
  
  
  // 约束10: 扩增子重叠约束（while循环合并版本）
  int overlap_constraint_count = 0;
  std::cout<<"collet constratins"<<std::endl;
  
  // 第一步：收集所有原始约束字符串，按D_var分组
  //unordered_map<string, vector<string>> constraints_by_dvar;
  unordered_map<string, unordered_map<string, unordered_set<string>>> pre_merged_constraints;
  
  for(int k = 0; k < num_loci; k++) {
    for(int l = k + 1; l < num_loci; l++) {
      string locus_k = loci[k];
      string locus_l = loci[l];
      
      std::cout<<"locus_k: "<<locus_k<<" locus_l: "<<locus_l<<std::endl;
      
      
      // 获取位点k和l的5'和3'引物
      if(locus_to_primers_5.find(locus_k) == locus_to_primers_5.end() ||
         locus_to_primers_3.find(locus_k) == locus_to_primers_3.end() ||
         locus_to_primers_5.find(locus_l) == locus_to_primers_5.end() ||
         locus_to_primers_3.find(locus_l) == locus_to_primers_3.end()) {
        continue;
      }
      
      const vector<int>& primers_5_k = locus_to_primers_5[locus_k];
      const vector<int>& primers_3_k = locus_to_primers_3[locus_k];
      const vector<int>& primers_5_l = locus_to_primers_5[locus_l];
      const vector<int>& primers_3_l = locus_to_primers_3[locus_l];
      
      // 遍历所有可能的引物组合
      for(int idx_5_k : primers_5_k) {
        for(int idx_3_k : primers_3_k) {
          for(int idx_5_l : primers_5_l) {
            for(int idx_3_l : primers_3_l) {
              // 检查这四个引物之间是否都兼容
              bool all_compatible = true;
              
              vector<pair<int, int>> cross_locus_pairs = {
                {idx_5_k, idx_5_l}, {idx_5_k, idx_3_l},
                {idx_3_k, idx_5_l}, {idx_3_k, idx_3_l}
              };
              
              for(auto& p : cross_locus_pairs) {
                int min_idx = min(p.first, p.second);
                int max_idx = max(p.first, p.second);
                if(compatible_set.find(make_pair(min_idx, max_idx)) == compatible_set.end()) {
                  all_compatible = false;
                  break;
                }
              }
              
              if(!all_compatible) continue;
              
              // 计算扩增子大小并检查重叠
              double amplicon_low_k = all_primers[idx_5_k].distance + 
                min(all_primers[idx_5_k].copy_min, all_primers[idx_3_k].copy_min) +
                all_primers[idx_3_k].distance;
              
              double amplicon_high_k = all_primers[idx_5_k].distance + 
                max(all_primers[idx_5_k].copy_max, all_primers[idx_3_k].copy_max) +
                all_primers[idx_3_k].distance;
              
              double amplicon_low_l = all_primers[idx_5_l].distance + 
                min(all_primers[idx_5_l].copy_min, all_primers[idx_3_l].copy_min) +
                all_primers[idx_3_l].distance;
              
              double amplicon_high_l = all_primers[idx_5_l].distance + 
                max(all_primers[idx_5_l].copy_max, all_primers[idx_3_l].copy_max) +
                all_primers[idx_3_l].distance;
              
              bool could_overlap = !(amplicon_high_k + size_gap <= amplicon_low_l ||
                                     amplicon_high_l + size_gap <= amplicon_low_k);
              
              if(!could_overlap) continue;
              
              // 创建约束字符串
              string d_var = "D" + to_string(k) + "_" + to_string(l);
              string x_5_k = "X" + to_string(idx_5_k);
              string x_3_k = "X" + to_string(idx_3_k);
              string x_5_l = "X" + to_string(idx_5_l);
              string x_3_l = "X" + to_string(idx_3_l);
              
              // 前三个X的组合作为键（排序以确保一致性）
              vector<string> key_parts = {x_5_k, x_3_k, x_5_l};
              sort(key_parts.begin(), key_parts.end());
              string key = key_parts[0] + "_" + key_parts[1] + "_" + key_parts[2];
              
              // 将第四个X添加到对应的集合中
              pre_merged_constraints[d_var][key].insert(x_3_l);
              
              
            }
          }
        }
      }
    }
  }
  
  std::cout<<"process constratins"<<std::endl;
  
  
  // 第二步：将预合并的结果转换为最终的约束字符串
  unordered_map<string, vector<string>> constraints_by_dvar;
  
  for(const auto& dvar_entry : pre_merged_constraints) {
    const string& d_var = dvar_entry.first;
    const auto& key_map = dvar_entry.second;
    
    for(const auto& key_entry : key_map) {
      const string& key = key_entry.first;
      const unordered_set<string>& x3_l_set = key_entry.second;
      
      // 解析键中的前三个X变量
      // 假设key格式: Xa_Xb_Xc
      vector<string> key_vars;
      size_t start = 0, end;
      while((end = key.find('_', start)) != string::npos) {
        key_vars.push_back(key.substr(start, end - start));
        start = end + 1;
      }
      key_vars.push_back(key.substr(start));
      
      if(key_vars.size() != 3) continue;
      
      // 创建合并后的约束
      if(x3_l_set.size() == 1) {
        // 只有一个第四个X，创建原始约束
        string constraint = d_var + " + " +
          key_vars[0] + " + " +
          key_vars[1] + " + " +
          key_vars[2] + " + " +
          *x3_l_set.begin() + " <= 4";
          
          constraints_by_dvar[d_var].push_back(constraint);
      } else {
        // 多个第四个X，创建合并约束
        string constraint = d_var + " + " +
          key_vars[0] + " + " +
          key_vars[1] + " + " +
          key_vars[2] + " + ";
        
        // 添加所有第四个X
        vector<string> sorted_x3_l(x3_l_set.begin(), x3_l_set.end());
        sort(sorted_x3_l.begin(), sorted_x3_l.end());
        
        for(size_t i = 0; i < sorted_x3_l.size(); i++) {
          if(i > 0) constraint += " + ";
          constraint += sorted_x3_l[i];
        }
        
        constraint += " <= " + to_string(3 + sorted_x3_l.size());
        
        constraints_by_dvar[d_var].push_back(constraint);
      }
    }
  }
  
  // 第二步：对每个D_var的约束进行合并
  // 在约束10开始之前，添加这些辅助函数作为lambda
  unordered_map<string, string> var_to_category;
  
  // 填充查找表（在约束10之前）
  for(int i = 0; i < num_primers; i++) {
    string var_name = "X" + to_string(i);
    var_to_category[var_name] = all_primers[i].locus + "_" + all_primers[i].end_type;
  }
  
  // 辅助函数：从约束字符串中提取 X 变量的整数索引集合
  auto extract_x_indices = [](const string& constraint) -> unordered_set<int> {
    unordered_set<int> indices;
    stringstream ss(constraint);
    string token;
    
    while(getline(ss, token, '+')) {
      // 清理空格和<=部分
      size_t pos = token.find('<');
      if(pos != string::npos) {
        token = token.substr(0, pos);
      }
      
      // 移除空格
      token.erase(remove_if(token.begin(), token.end(), ::isspace), token.end());
      
      if(!token.empty() && token[0] == 'X') {
        try {
          int idx = stoi(token.substr(1));
          indices.insert(idx);
        } catch(...) {
          // 忽略解析错误
        }
      }
      // 忽略 D 变量
    }
    return indices;
  };
  

  
  
  // 定义节点数据结构，每类变量分开存储
  // 定义节点数据结构，每类变量分开存储
  struct NodeCategories {
    bitset<MAX_PRIMERS> locus_k_5;
    bitset<MAX_PRIMERS> locus_k_3;
    bitset<MAX_PRIMERS> locus_l_5;
    bitset<MAX_PRIMERS> locus_l_3;
    int count_k_5;  // 预计算的大小
    int count_k_3;
    int count_l_5;
    int count_l_3;
    int first_k_5;  // 第一个引物索引（-1表示空）
    int first_k_3;
    int first_l_5;
    int first_l_3;
    string d_var;
    int k_index;
    int l_index;
    
    // 辅助函数：从bitset获取第一个设置位的索引
    int get_first_from_bitset(const bitset<MAX_PRIMERS>& bs) {
      if(bs.any()) {
        for(size_t i = 0; i < MAX_PRIMERS; ++i) {
          if(bs.test(i)) return static_cast<int>(i);
        }
      }
      return -1;
    }
    
    // 更新所有预计算值
    void update_all() {
      count_k_5 = (int)locus_k_5.count();
      count_k_3 = (int)locus_k_3.count();
      count_l_5 = (int)locus_l_5.count();
      count_l_3 = (int)locus_l_3.count();
      first_k_5 = get_first_from_bitset(locus_k_5);
      first_k_3 = get_first_from_bitset(locus_k_3);
      first_l_5 = get_first_from_bitset(locus_l_5);
      first_l_3 = get_first_from_bitset(locus_l_3);
    }
  };
  
  // 辅助函数：检查两个集合是否完全相同
  // 辅助函数：检查两个集合是否完全相同（改为lambda函数）
  auto sets_equal = [](const bitset<MAX_PRIMERS>& set1, 
                       const bitset<MAX_PRIMERS>& set2) -> bool {
                         return set1 == set2;  // bitset的比较是逐字（word）比较，非常快
                       };
  // 获取变量的类别（位点_端类型） - 使用整数索引版本
  auto get_var_category_simple = [&](int idx) -> string {
    if(idx < 0 || idx >= num_primers) return "";
    return all_primers[idx].locus + "_" + all_primers[idx].end_type;
  };
  
  // 辅助函数：从bitset获取第一个设置位的索引
  auto get_first_bit_idx = [&](const bitset<MAX_PRIMERS>& bs) -> int {
    if(bs.any()) {
      // 更高效的查找方法：使用_bittestandset等内在函数，或者简单循环
      for(size_t i = 0; i < num_primers; ++i) {  // 只检查实际引物数量
        if(bs.test(i)) return static_cast<int>(i);
      }
    }
    return -1;
  };
  
  auto get_edge_type_with_categories = [&](const NodeCategories& node1,
                                           const NodeCategories& node2) -> string {

                                             // 2. 使用预计算的计数检查大小差异
                                             int size_diff_count = 0;
                                             size_diff_count += (node1.count_k_5 != node2.count_k_5) ? 1 : 0;
                                             size_diff_count += (node1.count_k_3 != node2.count_k_3) ? 1 : 0;
                                             size_diff_count += (node1.count_l_5 != node2.count_l_5) ? 1 : 0;
                                             size_diff_count += (node1.count_l_3 != node2.count_l_3) ? 1 : 0;
                                             
                                             if (size_diff_count >= 2) {
                                               return "incompatible";
                                             }
                                             
                                             // 3. 检查具体哪一类不同
                                             int diff_count = 0;
                                             string diff_type = "";
                                             
                                             // 比较 locus_k_5
                                             if (node1.locus_k_5 != node2.locus_k_5) {
                                               diff_count++;
                                               if (diff_count > 1) return "incompatible";
                                               // 使用预存储的第一个引物索引
                                               int idx = (node1.first_k_5 != -1) ? node1.first_k_5 : node2.first_k_5;
                                               if (idx != -1) diff_type = get_var_category_simple(idx);
                                             }
                                             
                                             // 比较 locus_k_3
                                             if (node1.locus_k_3 != node2.locus_k_3) {
                                               diff_count++;
                                               if (diff_count > 1) return "incompatible";
                                               int idx = (node1.first_k_3 != -1) ? node1.first_k_3 : node2.first_k_3;
                                               if (idx != -1) diff_type = get_var_category_simple(idx);
                                             }
                                             
                                             // 比较 locus_l_5
                                             if (node1.locus_l_5 != node2.locus_l_5) {
                                               diff_count++;
                                               if (diff_count > 1) return "incompatible";
                                               int idx = (node1.first_l_5 != -1) ? node1.first_l_5 : node2.first_l_5;
                                               if (idx != -1) diff_type = get_var_category_simple(idx);
                                             }
                                             
                                             // 比较 locus_l_3
                                             if (node1.locus_l_3 != node2.locus_l_3) {
                                               diff_count++;
                                               if (diff_count > 1) return "incompatible";
                                               int idx = (node1.first_l_3 != -1) ? node1.first_l_3 : node2.first_l_3;
                                               if (idx != -1) diff_type = get_var_category_simple(idx);
                                             }
                                             
                                             // 1. 快速检查完全相同
                                             if (node1.locus_k_5 == node2.locus_k_5 &&
                                                 node1.locus_k_3 == node2.locus_k_3 &&
                                                 node1.locus_l_5 == node2.locus_l_5 &&
                                                 node1.locus_l_3 == node2.locus_l_3) {
                                               return "identical";
                                             }
                                             
                                             // 如果执行到这里，说明最多只有一类不同
                                             return diff_type;  // 可能为空（如果diff_count==0，但理论上不应该）
                                           };
  
  auto extract_vars = [](const string& constraint) -> unordered_set<string> {
    unordered_set<string> vars;
    stringstream ss(constraint);
    string token;
    
    while(getline(ss, token, '+')) {
      // 清理空格和<=部分
      size_t pos = token.find('<');
      if(pos != string::npos) {
        token = token.substr(0, pos);
      }
      
      // 移除空格
      token.erase(remove_if(token.begin(), token.end(), ::isspace), token.end());
      
      if(!token.empty() && (token[0] == 'D' || token[0] == 'X')) {
        vars.insert(token);
      }
    }
    return vars;
  };
  // 格式化约束字符串（添加 D 变量）
  auto format_constraint_with_dvar = [&](const unordered_set<int>& x_indices, 
                                         const string& d_var) -> string {
                                           // 创建变量名列表
                                           vector<string> all_vars;
                                           all_vars.push_back(d_var);
                                           
                                           for(int idx : x_indices) {
                                             all_vars.push_back("X" + to_string(idx));
                                           }
                                           
                                           // 排序（D变量在前，X变量按数字排序）
                                           sort(all_vars.begin() + 1, all_vars.end()); // 只排序X变量
                                           
                                           // 构建约束字符串
                                           string constraint = "";
                                           for(size_t i = 0; i < all_vars.size(); i++) {
                                             if(i > 0) constraint += " + ";
                                             constraint += all_vars[i];
                                           }
                                           constraint += " <= " + to_string(all_vars.size() - 1);
                                           
                                           return constraint;
                                         };
  
  vector<string> final_constraints;
  
  int locus_index=1;
  for(const auto& dvar_entry : constraints_by_dvar) {
    const string& d_var = dvar_entry.first;
    const vector<string>& d_constraints = dvar_entry.second;
    
    std::cout << "d_var: " << d_var << " "<<locus_index++<< std::endl;
    
    if(d_constraints.empty()) continue;
    
    // 1. 初始化：每个约束作为一个节点（只包含X变量索引）
    // 初始化：每个约束转换为NodeCategories结构
    vector<NodeCategories> categorized_nodes;
    categorized_nodes.reserve(d_constraints.size());
    
    for(const string& c : d_constraints) {
      unordered_set<int> x_indices = extract_x_indices(c);
      
      NodeCategories node;
      node.d_var = d_var;
      
      // 解析 d_var 获取位点索引
      string d_num = d_var.substr(1);
      size_t underscore_pos = d_num.find('_');
      if(underscore_pos == string::npos) continue;
      
      node.k_index = stoi(d_num.substr(0, underscore_pos));
      node.l_index = stoi(d_num.substr(underscore_pos + 1));
      
      string locus_k = loci[node.k_index];
      string locus_l = loci[node.l_index];
      
      for(int idx : x_indices) {
        const PrimerInfo& primer = all_primers[idx];
        
        if(primer.locus == locus_k) {
          if(primer.end_type == "5") {
            node.locus_k_5.set(idx);  // 使用 set() 而不是 insert()
          } else if(primer.end_type == "3") {
            node.locus_k_3.set(idx);
          }
        } else if(primer.locus == locus_l) {
          if(primer.end_type == "5") {
            node.locus_l_5.set(idx);
          } else if(primer.end_type == "3") {
            node.locus_l_3.set(idx);
          }
        }
      }
      
      node.update_all();
      
      categorized_nodes.push_back(node);
    }
    
    
    int max_rounds = 5;
    
    // 2. 多轮合并
    for(int round = 0; round < max_rounds; round++) {
      int n = categorized_nodes.size();
      
      std::cout << "round: " << round << " n= " << n ;
      
      if(n <= 1) break;  // 没有可以合并的了
      
      // 2.1 计算所有节点对之间的边类型
      //vector<vector<string>> edge_types(n, vector<string>(n, ""));
      //vector<vector<bool>> can_merge(n, vector<bool>(n, false));
      
      std::cout << " get edge types "  ;
      if(n > 0) {
        std::cout << " node size: " << categorized_nodes[0].locus_l_3.size() << std::endl ;
      } else {
        std::cout << std::endl;
      }
      
      vector<int> edge_count(n, 0);  // 每个节点的边数
      vector<unordered_set<string>> node_edge_types(n);  // 每个节点已有的边类型
      vector<vector<pair<int, string>>> edges(n);  // 每个节点的边列表 (邻居, 边类型)
      
      int edge_type_count_thres=1;
      
      
      
      // 只计算必要的边
      for(int i = 0; i < n; i++) {
        if(edge_count[i] >= edge_type_count_thres) continue;
        
        NodeCategories& node_i = categorized_nodes[i];
        
        for(int j = i + 1; j < n; j++) {
          if(edge_count[j] >= edge_type_count_thres) continue;
          
          NodeCategories& node_j = categorized_nodes[j];
          
          // 直接使用新的get_edge_type_with_categories
          string type = get_edge_type_with_categories(node_i, node_j);
          //string type="xxx";
          // 检查是否可以合并
          if(type == "incompatible" || type.empty()) continue;
          
          // 检查双方是否已经有这种边类型
          if(node_edge_types[i].find(type) != node_edge_types[i].end()) continue;
          if(node_edge_types[j].find(type) != node_edge_types[j].end()) continue;
          
          // 添加边
          edges[i].push_back(make_pair(j, type));
          edges[j].push_back(make_pair(i, type));
          
          node_edge_types[i].insert(type);
          node_edge_types[j].insert(type);
          
          edge_count[i]++;
          edge_count[j]++;
          
          if(edge_count[i] >= edge_type_count_thres) break;
        }
      }
      std::cout << "get edge types end" << std::endl;
      
      // 2.2 按边类型分别进行合并
      vector<bool> merged(n, false);
      //vector<unordered_set<int>> next_round_nodes;
      vector<NodeCategories> next_round_nodes;
      
      // 收集所有存在的边类型
      unordered_set<string> all_edge_types;
      for(int i = 0; i < n; i++) {
        for(const auto& edge : edges[i]) {
          all_edge_types.insert(edge.second);
        }
      }
      
      // 对每种边类型分别处理
      for(const string& edge_type : all_edge_types) {
        // 为当前边类型创建连通分量
        for(int i = 0; i < n; i++) {
          if(merged[i]) continue;
          
          // 检查节点i是否有当前类型的边
          bool has_this_type = false;
          for(const auto& edge : edges[i]) {
            if(edge.second == edge_type) {
              has_this_type = true;
              break;
            }
          }
          
          if(!has_this_type) continue;
          
          // 收集当前连通分量
          vector<int> component;
          unordered_set<int> visited;
          queue<int> q;
          
          q.push(i);
          visited.insert(i);
          
          while(!q.empty()) {
            int curr = q.front();
            q.pop();
            component.push_back(curr);
            
            // 只遍历当前边类型的邻居
            for(const auto& edge : edges[curr]) {
              int neighbor = edge.first;
              string type = edge.second;
              
              if(type == edge_type && 
                 visited.find(neighbor) == visited.end() && 
                 !merged[neighbor]) {
                 // 检查与组件内所有节点的兼容性
                 bool compatible_with_all = true;
                for(int existing : component) {
                  bool found = false;
                  for(const auto& e : edges[existing]) {
                    if(e.first == neighbor && e.second == edge_type) {
                      found = true;
                      break;
                    }
                  }
                  if(!found) {
                    compatible_with_all = false;
                    break;
                  }
                }
                
                if(compatible_with_all) {
                  q.push(neighbor);
                  visited.insert(neighbor);
                }
              }
            }
          }
          
          // 标记为已合并
          for(int node_idx : component) {
            merged[node_idx] = true;
          }
          
          // 合并当前连通分量
          if(component.size() > 1) {
            NodeCategories merged_node;
            merged_node.d_var = categorized_nodes[component[0]].d_var;
            merged_node.k_index = categorized_nodes[component[0]].k_index;
            merged_node.l_index = categorized_nodes[component[0]].l_index;
            
            for(int node_idx : component) {
              const NodeCategories& node = categorized_nodes[node_idx];
              
              // 合并 locus_k_5（按位或）
              merged_node.locus_k_5 |= node.locus_k_5;
              // 合并 locus_k_3
              merged_node.locus_k_3 |= node.locus_k_3;
              // 合并 locus_l_5
              merged_node.locus_l_5 |= node.locus_l_5;
              // 合并 locus_l_3
              merged_node.locus_l_3 |= node.locus_l_3;
            }
            
            // 关键：重新计算合并后的计数
            merged_node.update_all();
            
            
            next_round_nodes.push_back(merged_node);
            
          } else {
            // 只有一个节点，取消合并标记（可能被其他类型合并）
            merged[component[0]] = false;
          }
        }
      }
      
      // 处理剩余的孤立节点（没有边或边类型为incompatible）
      for(int i = 0; i < n; i++) {
        if(!merged[i]) {
          next_round_nodes.push_back(categorized_nodes[i]);
        }
      }
      
      // 重置edge_count用于下一轮
      for(int i = 0; i < n; i++) {
        edge_count[i] = 0;
      }
      
      // 准备下一轮
      categorized_nodes = next_round_nodes;
      
      Rcpp::Rcout << "  Round " << (round + 1) << ": " 
                  << n << " -> " << categorized_nodes.size() 
                  << " constraints" << endl;
      
      if(categorized_nodes.size() <= 1) {
        break;  // 完全合并了
      }
    }
    
    // 3. 输出最终约束（添加回D变量）
    for(const auto& node : categorized_nodes) {
      // 合并四类变量为一个集合
      bitset<MAX_PRIMERS> all_indices;
      all_indices = node.locus_k_5 | node.locus_k_3 | node.locus_l_5 | node.locus_l_3;
      
      // 将 bitset 转换为整数索引集合
      unordered_set<int> indices_set;
      for(size_t i = 0; i < num_primers; ++i) {
        if(all_indices.test(i)) {
          indices_set.insert(static_cast<int>(i));
        }
      }
      
      string constraint = format_constraint_with_dvar(indices_set, node.d_var);
      final_constraints.push_back(constraint);
    }
    
  }
  
  // 第三步：去重并输出最终约束
  overlap_constraint_count = 0;
  unordered_set<string> unique_constraint_set;
  
  for(const string& c : final_constraints) {
    // 提取变量并排序，使约束标准化
    unordered_set<string> vars = extract_vars(c);
    vector<string> sorted_vars(vars.begin(), vars.end());
    sort(sorted_vars.begin(), sorted_vars.end());
    
    // 创建标准化字符串
    string normalized = "";
    for(size_t i = 0; i < sorted_vars.size(); i++) {
      if(i > 0) normalized += " + ";
      normalized += sorted_vars[i];
    }
    normalized += " <= 4";
    
    // 检查是否已存在
    if(unique_constraint_set.find(normalized) == unique_constraint_set.end()) {
      unique_constraint_set.insert(normalized);
      
      // 输出约束
      lp << normalized << "\n";
      constr_count++;
      overlap_constraint_count++;
    }
  }
  // 输出统计信息
  int total_original = 0;
  for(const auto& entry : constraints_by_dvar) {
    total_original += entry.second.size();
  }
  
  Rcpp::Rcout << "Overlap constraints summary:" << endl;
  Rcpp::Rcout << "  Original constraints: " << total_original << endl;
  Rcpp::Rcout << "  Final constraints: " << overlap_constraint_count << endl;
  if(total_original > 0) {
    double reduction = (1.0 - (double)overlap_constraint_count / total_original) * 100;
    Rcpp::Rcout << "  Reduction: " << fixed << setprecision(1) << reduction << "%" << endl;
  }
  
  
  // 变量类型定义
  lp << "\nBinary\n";
  
  // X变量（引物选择）
  for(int i = 0; i < num_primers; i++) {
    lp << "X" << i << "\n";
  }
  
  // P变量（位点选中）
  for(const string& p_var : P_vars) {
    lp << p_var << "\n";
  }
  
  // C变量（位点-通道分配）
  for(const string& c_var : C_vars) {
    lp << c_var << "\n";
  }
  
  // D变量（位点之间的）
  for(const string& d_var : D_vars) {
    lp << d_var << "\n";
  }
  
  // E变量（辅助变量）
  for(const string& e_var : E_vars) {
    lp << e_var << "\n";
  }
  
  // Z变量
  for(const string& z_var : Z_vars) {
    lp << z_var << "\n";
  }
  
  lp << "\nEnd\n";
  lp.close();
  
  // 4. 保存引物信息
  string primers_filename = "single_primers_info.tsv";
  ofstream primers_file(primers_filename);
  primers_file << "primer_index\tprimer_id\tlocus\tlocus_index\tend_type\tprimer_seq\tdistance\tcopy_min\tcopy_max\trepeat_distance\tobjective_coefficient" << endl;
  
  for(int i = 0; i < num_primers; i++) {
    double coefficient = 50000.0 - all_primers[i].distance * 1.0;
    
    primers_file << i << "\t"
                 << all_primers[i].primer_id << "\t"
                 << all_primers[i].locus << "\t"
                 << locus_to_index[all_primers[i].locus] << "\t"
                 << all_primers[i].end_type << "\t"
                 << all_primers[i].primer_seq << "\t"
                 << all_primers[i].distance << "\t"
                 << all_primers[i].copy_min << "\t"
                 << all_primers[i].copy_max << "\t"
                 << all_primers[i].repeat_dist << "\t"
                 << coefficient << endl;
  }
  primers_file.close();
  
  // 预览LP文件
  ifstream preview_file(lp_filename);
  string lp_preview(1000, ' ');
  preview_file.read(&lp_preview[0], 1000);
  preview_file.close();
  
  // 统计信息
  int total_incompatible_pairs = 0;
  for(int i = 0; i < num_primers; i++) {
    for(int j = i + 1; j < num_primers; j++) {
      if(all_primers[i].locus == all_primers[j].locus) continue;
      
      bool is_compatible = (compatible_set.find(make_pair(i, j)) != compatible_set.end());
      if(!is_compatible) {
        total_incompatible_pairs++;
      }
    }
  }
  
  // 返回结果
  List result = List::create(
    Named("primers_info_file") = primers_filename,
    Named("lp_file") = lp_filename,
    Named("num_variables") = X_vars.size() + P_vars.size() + C_vars.size() + 
      D_vars.size() + E_vars.size() + Z_vars.size(),
      Named("num_constraints") = constr_count,
      Named("num_primers") = num_primers,
      Named("num_loci") = num_loci,
      Named("num_channels") = num_channels,
      Named("compatible_pairs") = (int)compatible_set.size(),
      Named("incompatible_pairs") = total_incompatible_pairs,
      Named("overlap_constraints") = overlap_constraint_count,
      Named("model_summary") = List::create(
        Named("X_vars") = (int)X_vars.size(),
        Named("P_vars") = (int)P_vars.size(),
        Named("C_vars") = (int)C_vars.size(),
        Named("D_vars") = (int)D_vars.size(),
        Named("E_vars") = (int)E_vars.size(),
        Named("Z_vars") = (int)Z_vars.size()
      ),
      Named("lp_preview") = lp_preview + "..."
  );
  
  return result;
}