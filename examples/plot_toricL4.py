import matplotlib.pyplot as plt
import numpy as np

# 数据
physical_error_rate = [0.005, 0.01, 0.015, 0.02, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4]
logical_error_rate_ML = [0.0361328, 0.0673828, 0.108398, 0.277344, 0.274414, 0.297852, 0.537109, 0.707031, 0.75293, 0.734375]
logical_error_rate_mode0 = [0, 0.00292969, 0.0078125, 0.0126953, 0.0175781, 0.0761719, 0.267578, 0.614258, 0.756836, 0.756836]

# 创建图形和坐标轴
plt.figure(figsize=(10, 6))

# 绘制双对数图
plt.plot(physical_error_rate, logical_error_rate_ML, '-^', markersize=8, linewidth=2, label='ML')
plt.plot(physical_error_rate, logical_error_rate_mode0, '-o', markersize=8, linewidth=2, label='Mode 0')

# 设置标题和标签
plt.title('Logical Error Rate vs Physical Error Rate', fontsize=16)
plt.xlabel('Physical Error Rate', fontsize=14)
plt.ylabel('Logical Error Rate', fontsize=14)

# 设置网格
plt.grid(True, which="both", ls="-", alpha=0.2)

# 调整坐标轴
plt.xlim(0.004, 0.5)
plt.ylim(0.0001, 1)

# 显示图例
plt.legend(loc='best')

# 保存图形
plt.savefig('error_rate_plot.png', dpi=300, bbox_inches='tight')

# 显示图形
plt.show()