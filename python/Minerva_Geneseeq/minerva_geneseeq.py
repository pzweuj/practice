#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MINERVA评分计算程序
https://minerva.geneseeq.com/
https://pubmed.ncbi.nlm.nih.gov/34750392/
MINERVA评分用于评估EGFR突变阳性非小细胞肺癌患者的辅助治疗策略
"""

from typing import List, Dict, Tuple
from dataclasses import dataclass
from enum import Enum
import argparse
import sys


class MinervaGroup(Enum):
    """MINERVA评分分组"""
    HTP = 0  # Highly TKI-Preferable group：强烈推荐靶向治疗
    TP = 1   # TKI-Preferable group：推荐靶向治疗
    CP = 2   # Chemo-Preferable Group：推荐化疗


@dataclass
class MinervaInput:
    """MINERVA评分输入参数"""
    cancer_staging: str = ""       # 癌症分期
    egfr: str = "not_detected"     # EGFR突变状态 (19del/L858R/other/not_detected)
    rb1: bool = False              # RB1基因状态
    nkx2: bool = False             # NKX2-1基因状态
    tp53: bool = False             # TP53基因状态
    myc: bool = False              # MYC基因状态
    cdk4: bool = False             # CDK4基因状态
    
    def __post_init__(self):
        """初始化默认值"""
        pass


@dataclass
class MinervaResult:
    """MINERVA评分结果"""
    score: float                         # 总评分
    group: MinervaGroup                  # 分组
    group_title: str                     # 分组标题
    group_description: str               # 分组描述
    individual_scores: Dict[str, float]  # 各基因的评分


class MinervaCalculator:
    """MINERVA评分计算器"""
    
    # 基因评分权重（来自html中的calc方法）
    GENE_SCORES = {
        'rb1': 2.88,    # RB1阳性时的评分
        'nkx2': -2.72,  # NKX2-1阳性时的评分
        'cdk4': -2.26,  # CDK4阳性时的评分
        'tp53': -2.11,  # TP53阳性时的评分
        'myc': -1.98,   # MYC阳性时的评分
    }
    
    # 分组阈值
    GROUP_THRESHOLDS = {
        'htp_upper': -0.5,  # HTP组上限
        'tp_upper': 0.5,    # TP组上限
    }
    
    # 分组信息
    GROUP_INFO = {
        MinervaGroup.HTP: {
            'title': 'HTP组（Highly TKI-Preferable group）：强烈推荐靶向治疗',
            'description': 'HTP组人群更适合使用EGFR-TKI作为辅助治疗（HR=0.21）：吉非替尼组中位无病生存期（mDFS）为34.5个月，显著优于化疗组9.1个月；5年生存率辅助吉非替尼治疗约67.3%，对比辅助化疗约38.3%。'
        },
        MinervaGroup.TP: {
            'title': 'TP组（TKI-Preferable group）：推荐靶向治疗',
            'description': 'TP组人群使用吉非替尼辅助治疗的获益情况与ADJUVANT/CTONG1104临床试验的意向治疗人群相似 (HR=0.61)：吉非替尼组中位无病生存期（mDFS）为32.8个月，优于化疗组20.7个月。'
        },
        MinervaGroup.CP: {
            'title': 'CP组（Chemo-Preferable Group）：推荐化疗',
            'description': 'CP组人群在携带EGFR敏感突变的情况下，使用EGFR-TKI进行辅助治疗的疗效可能仍劣于化疗（HR=3.06）。化疗组中位无病生存期（mDFS）为34.2个月，显著高于吉非替尼组19.3个月。5年生存率辅助化疗约61.5%，对比辅助吉非替尼治疗约28.3%。'
        }
    }
    
    def calculate(self, input_data: MinervaInput) -> MinervaResult:
        """
        计算MINERVA评分
        
        Args:
            input_data: 输入参数
            
        Returns:
            MinervaResult: 评分结果
        """
        # 验证必填项
        if not input_data.cancer_staging:
            raise ValueError("癌症分期为必填项")
        if not input_data.egfr or input_data.egfr == "not_detected":
            # EGFR未检出时，MINERVA评分不适用
            raise ValueError("EGFR突变状态为未检出，MINERVA评分仅适用于EGFR突变阳性患者")
        
        # 计算各基因评分
        individual_scores = {}
        
        # RB1评分：如果有RB1拷贝数缺失或突变，得分2.88，否则0
        rb1_score = self.GENE_SCORES['rb1'] if input_data.rb1 else 0.0
        individual_scores['rb1'] = rb1_score
        
        # NKX2评分：如果有NKX2突变，得分-2.72，否则0
        nkx2_score = self.GENE_SCORES['nkx2'] if input_data.nkx2 else 0.0
        individual_scores['nkx2'] = nkx2_score
        
        # CDK4评分：如果有CDK4突变，得分-2.26，否则0
        cdk4_score = self.GENE_SCORES['cdk4'] if input_data.cdk4 else 0.0
        individual_scores['cdk4'] = cdk4_score
        
        # TP53评分：如果有TP53突变，得分-2.11，否则0
        tp53_score = self.GENE_SCORES['tp53'] if input_data.tp53 else 0.0
        individual_scores['tp53'] = tp53_score
        
        # MYC评分：如果有MYC突变，得分-1.98，否则0
        myc_score = self.GENE_SCORES['myc'] if input_data.myc else 0.0
        individual_scores['myc'] = myc_score
        
        # 计算总评分
        total_score = rb1_score + nkx2_score + cdk4_score + tp53_score + myc_score
        
        # 根据评分确定分组
        group = self._determine_group(total_score)
        
        # 获取分组信息
        group_info = self.GROUP_INFO[group]
        
        return MinervaResult(
            score=round(total_score, 2),
            group=group,
            group_title=group_info['title'],
            group_description=group_info['description'],
            individual_scores=individual_scores
        )
    
    def _determine_group(self, score: float) -> MinervaGroup:
        """
        根据评分确定分组
        
        Args:
            score: 总评分
            
        Returns:
            MinervaGroup: 分组结果
        """
        if score <= self.GROUP_THRESHOLDS['htp_upper']:
            return MinervaGroup.HTP
        elif score < self.GROUP_THRESHOLDS['tp_upper']:
            return MinervaGroup.TP
        else:
            return MinervaGroup.CP


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="MINERVA评分计算器 - 用于评估EGFR突变阳性非小细胞肺癌患者的辅助治疗策略",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法（仅必填参数）
  python minerva_geneseeq.py --cancer-staging IIIA --egfr 19del
  
  # 完整用法（包含基因突变信息）
  python minerva_geneseeq.py --cancer-staging IB --egfr L858R --rb1 --tp53 --nkx2 --myc --cdk4
  
  # RB1缺失示例
  python minerva_geneseeq.py --cancer-staging IIIA --egfr 19del --rb1 --tp53
  
  # 其他EGFR突变示例
  python minerva_geneseeq.py --cancer-staging IIA --egfr other --myc --cdk4
  
  # 运行内置示例
  python minerva_geneseeq.py --demo
        """
    )
    
    # 必填参数组
    required = parser.add_argument_group('必填参数')
    required.add_argument(
        '--cancer-staging', 
        type=str, 
        choices=['IB', 'IIA', 'IIB', 'IIIA', 'IIIB', 'IV'],
        help='癌症分期: IB, IIA, IIB, IIIA, IIIB, IV (MINERVA评分适用于IB-IV期)'
    )
    required.add_argument(
        '--egfr', 
        type=str, 
        choices=['19del', 'L858R', 'other', 'not_detected'],
        default='not_detected',
        help='EGFR突变状态: 19del=19外显子缺失, L858R=L858R突变, other=其他突变, not_detected=未检出 (默认: not_detected)'
    )
    
    # 可选参数组 - 基因突变信息
    optional = parser.add_argument_group('可选参数 - 基因突变信息')
    optional.add_argument(
        '--rb1', 
        action='store_true',
        help='RB1基因阳性 (拷贝数缺失或突变)'
    )
    optional.add_argument(
        '--nkx2', 
        action='store_true',
        help='NKX2-1基因阳性 (存在拷贝数扩增)'
    )
    optional.add_argument(
        '--tp53', 
        action='store_true',
        help='TP53基因阳性 (NM_000546第四或第五外显子错义突变)'
    )
    optional.add_argument(
        '--myc', 
        action='store_true',
        help='MYC基因阳性 (存在拷贝数扩增)'
    )
    optional.add_argument(
        '--cdk4', 
        action='store_true',
        help='CDK4基因阳性 (存在拷贝数扩增)'
    )
    
    # 其他可选参数
    other = parser.add_argument_group('其他可选参数')
    other.add_argument(
        '--demo', 
        action='store_true',
        help='运行内置示例'
    )
    other.add_argument(
        '--verbose', 
        action='store_true',
        help='显示详细信息'
    )
    
    return parser.parse_args()


def run_demo():
    """运行内置示例"""
    print("=== MINERVA评分计算器 - 内置示例 ===\n")
    
    # 创建计算器实例
    calculator = MinervaCalculator()
    
    # 示例1：HTP组患者
    print("示例1：HTP组患者")
    input1 = MinervaInput(
        cancer_staging="IIIA",
        egfr="19del",
        rb1=True,     # RB1阳性，+2.88分
        nkx2=False,   # NKX2阴性，0分
        tp53=True,    # TP53阳性，-2.11分
        myc=False,    # MYC阴性，0分
        cdk4=False    # CDK4阴性，0分
    )
    
    result1 = calculator.calculate(input1)
    print(f"总评分: {result1.score}")
    print(f"分组: {result1.group_title}")
    print(f"建议: {result1.group_description}")
    print(f"各基因评分: {result1.individual_scores}")
    print()
    
    # 示例2：CP组患者
    print("示例2：CP组患者")
    input2 = MinervaInput(
        cancer_staging="IB",
        egfr="L858R",
        rb1=False,        # RB1阴性，0分
        nkx2=True,        # NKX2阳性，-2.72分
        tp53=True,        # TP53阳性，-2.11分
        myc=True,         # MYC阳性，-1.98分
        cdk4=True         # CDK4阳性，-2.26分
    )
    
    result2 = calculator.calculate(input2)
    print(f"总评分: {result2.score}")
    print(f"分组: {result2.group_title}")
    print(f"建议: {result2.group_description}")
    print(f"各基因评分: {result2.individual_scores}")
    print()
    
    # 示例3：TP组患者
    print("示例3：TP组患者")
    input3 = MinervaInput(
        cancer_staging="IIA",
        egfr="19del",
        rb1=False,    # 所有基因都阴性，总分为0，属于TP组
        nkx2=False,
        tp53=False,
        myc=False,
        cdk4=False
    )
    
    result3 = calculator.calculate(input3)
    print(f"总评分: {result3.score}")
    print(f"分组: {result3.group_title}")
    print(f"建议: {result3.group_description}")
    print(f"各基因评分: {result3.individual_scores}")


def main():
    """主函数"""
    # 如果没有提供任何参数，显示帮助信息
    if len(sys.argv) == 1:
        parser = argparse.ArgumentParser(
            description="MINERVA评分计算器 - 用于评估EGFR突变阳性非小细胞肺癌患者的辅助治疗策略",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
使用示例:
  # 基本用法（仅必填参数）
  python minerva_geneseeq.py --cancer-staging IIIA --egfr 19del
  
  # 完整用法（包含基因突变信息）
  python minerva_geneseeq.py --cancer-staging IB --egfr L858R --rb1 --tp53 --nkx2 --myc --cdk4
  
  # RB1缺失示例
  python minerva_geneseeq.py --cancer-staging IIIA --egfr 19del --rb1 --tp53
  
  # 其他EGFR突变示例
  python minerva_geneseeq.py --cancer-staging IIA --egfr other --myc --cdk4
  
  # 运行内置示例
  python minerva_geneseeq.py --demo
        """
        )
        parser.print_help()
        return
    
    args = parse_arguments()
    
    # 如果指定了demo参数，运行内置示例
    if args.demo:
        run_demo()
        return
    
    # 检查必填参数
    if not args.cancer_staging:
        print("错误: 缺少必填参数")
        print("请使用 --help 查看使用说明")
        print("\n最简使用示例:")
        print("python minerva_geneseeq.py --cancer-staging IIIA --egfr 19del")
        sys.exit(1)
    
    # 创建输入数据
    input_data = MinervaInput(
        cancer_staging=args.cancer_staging,
        egfr=args.egfr,
        rb1=args.rb1,
        nkx2=args.nkx2,
        tp53=args.tp53,
        myc=args.myc,
        cdk4=args.cdk4
    )
    
    # 创建计算器并计算结果
    calculator = MinervaCalculator()
    
    try:
        result = calculator.calculate(input_data)
        
        # 输出结果
        print("=== MINERVA评分结果 ===")
        print(f"癌症分期: {input_data.cancer_staging}")
        egfr_status = {'19del': '19外显子缺失', 'L858R': 'L858R突变', 'other': '其他突变', 'not_detected': '未检出'}.get(input_data.egfr, input_data.egfr)
        print(f"EGFR状态: {egfr_status}")
        
        if args.verbose:
            print("\n基因突变状态:")
            print(f"  RB1: {'阳性' if input_data.rb1 else '阴性'}")
            print(f"  NKX2-1: {'阳性' if input_data.nkx2 else '阴性'}")
            print(f"  TP53: {'阳性' if input_data.tp53 else '阴性'}")
            print(f"  MYC: {'阳性' if input_data.myc else '阴性'}")
            print(f"  CDK4: {'阳性' if input_data.cdk4 else '阴性'}")
            
            print(f"\n各基因评分: {result.individual_scores}")
        
        print(f"\n总评分: {result.score}")
        print(f"分组: {result.group_title}")
        print(f"\n治疗建议:")
        print(f"{result.group_description}")
        
    except ValueError as e:
        print(f"错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
