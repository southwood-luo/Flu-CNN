import matplotlib.pyplot as plt
import matplotlib.font_manager as fm  # 新增引入该库
import math
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

# 设置字体为微软雅黑
myfont = fm.FontProperties(size=5)  # 修改字体和大小


def label(xy, text, font_size=6, angle=0):
    y = xy[1] - 0.05
    plt.text(xy[0], y, text, ha="center", family='Arial', size=font_size, rotation=angle)


def draw_wedge(data, center, radius, start_angle, end_angle, model, gene, colorset='red'):
    n = len(data)
    x = []
    y = []
    for i in range(n):
        angle = (start_angle + (i + 2) * (end_angle - start_angle) / (n + 3)) / 180
        xi = center[0] + (radius + data[i]) * math.cos(math.pi * angle)
        yi = center[1] + (radius + data[i]) * math.sin(math.pi * angle)
        x.append(xi)
        y.append(yi)
    if gene == 'PB2':
        plt.scatter(x, y, color=colorset, s=3, alpha=0.8, label=model)
    else:
        plt.scatter(x, y, color=colorset, s=3, alpha=0.8)


def plot_line(rmax, center, radius, end_angle, visibility_flag='none'):
    angle = end_angle / 180
    x1 = center[0] + radius * math.cos(math.pi * angle)
    y1 = center[1] + radius * math.sin(math.pi * angle)
    x2 = center[0] + (rmax + radius) * math.cos(math.pi * angle)
    y2 = center[1] + (rmax + radius) * math.sin(math.pi * angle)
    if visibility_flag == 'none':
        plt.plot([x1, x2], [y1, y2], color='grey', linewidth=0.5, alpha=0.2)
    else:
        plt.plot([x1, x2], [y1, y2], linewidth=0.5, alpha=0.3)
    step = 0.1
    standard = 0
    while standard < rmax:
        if visibility_flag == 'none':
            break
        if 0.4 < angle < 0.9:
            x_value = center[0] + (standard + radius) * math.cos(math.pi * (angle + 0.03))
            y_value = center[1] + (standard + radius) * math.sin(math.pi * (angle + 0.03))
        elif angle <= 0.4:
            x_value = center[0] + (standard + radius) * math.cos(math.pi * (angle + 0.02))
            y_value = center[1] + (standard + radius) * math.sin(math.pi * (angle + 0.02))
        elif 0.9 <= angle <= 1.3:
            x_value = center[0] + (standard + radius) * math.cos(math.pi * (angle + 0.017))
            y_value = center[1] + (standard + radius) * math.sin(math.pi * (angle + 0.017))
        else:
            x_value = center[0] + (standard + radius) * math.cos(math.pi * angle)
            y_value = center[1] + (standard + radius) * math.sin(math.pi * angle)
        plt.text(round(x_value, 5), round(y_value, 5), str(standard), fontproperties=myfont)  # 修改字体为微软雅黑
        standard = round(standard + step, 2)


def swedges(ary, center, r, pp):
    n = len(ary)
    d = float(360.0 / n)
    xd = d * 0.15
    st = 0
    en = 0 + d * 0.85
    ppdict = dict()
    ppw = []
    for a in ary:
        if pp[a - 1] not in ppdict:
            ppdict[pp[a - 1]] = [center, r, st, en]
        wedge = mpatches.Wedge(center, r, st, en, 0.12, ec="none")
        labelx = center[0] + 0.35 * math.cos(math.pi * (((st + en) / 2) / 180))
        labely = center[1] + 0.025 + 0.35 * math.sin(math.pi * (((st + en) / 2) / 180))
        label([labelx, labely], pp[a - 1])
        st = en + xd
        en = st + d * 0.85
        ppw.append(wedge)
    return ppw, ppdict


def loaddata(datacsv):
    f1 = open(datacsv)
    lines1 = f1.readlines()
    f1.close()
    datadict = dict()
    for line in lines1:
        if 'Sample Name' in line:
            continue
        if ',' in line:
            infs = line.strip().split(',')
            if len(infs) < 5:
                continue
            geneinf = infs[2]
            modelinf = infs[3]
            rate = float(infs[4])
            if geneinf in datadict:
                if modelinf in datadict[geneinf]:
                    datadict[geneinf][modelinf].append(rate)
                else:
                    datadict[geneinf][modelinf] = [rate]
            else:
                datadict[geneinf] = dict()
                datadict[geneinf][modelinf] = [rate]
    return datadict


subtypes = ['h9n2', 'h5n1', 'h7n9']
# color_set = {'VIDHOP': '#1f77b4', 'ML-DNTs': '#ff7f0e', 'Flu-CNN': '#d62728', 'FluPhenotype': '#2ca02c'}  # 修改颜色
color_set = {'VIDHOP': '#1f77b4', 'ML-DNTs': '#ff7f0e', 'Flu-CNN': '#d62728', 'FluPhenotype': '#2ca02c','Phylogenic': '#9467bd'}  # 修改颜色
for subtype in subtypes:
    grid = np.mgrid[0.2:0.8:3j, 0.2:0.8:3j].reshape(2, -1).T
    fig, ax = plt.subplots()
    patches = []
    numbers = range(1, 9)
    pp = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
    pp, ppdict = swedges(numbers, grid[2], 0.4, pp)
    patches.extend(pp)

    label(grid[2], subtype.upper(), 16)

    data_csv = f'hostdetection_res/hostdetection_mothods_re_{subtype}2.csv'
    data_dict = loaddata(data_csv)
    nature_colors = ['#74add1', '#fc8d62', '#8da0cb', '#e78ac3',
                     '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3']
    wedges = []
    RF = 0
    for gene in ppdict:
        center, radius, start_angle, end_angle = ppdict[gene]
        if gene not in data_dict:
            print('Error in data')
        else:
            n_model = len(data_dict[gene])
            model_list = list(data_dict[gene].keys())
            for i in range(n_model):
                rate_list = data_dict[gene][model_list[i]]
                ai = max(rate_list)
                width = (end_angle - start_angle) / n_model
                start_i = start_angle + i * width
                end_i = start_i + width
                wedge = mpatches.Wedge(center, radius + ai, start_i, end_i, ai, ec="none")
                RF = RF + 1
                wedges.append(wedge)
                draw_wedge(rate_list, center, radius, start_i, end_i, model_list[i], gene, color_set[model_list[i]])
                if model_list[i] == 'Flu-CNN':
                    plot_line(ai, center, radius, end_i, 'show')
                else:
                    plot_line(ai, center, radius, end_i, 'none')

    collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=0.8)
    colors = np.array(nature_colors)
    collection.set_facecolors(colors)
    ax.add_collection(collection)

    collection1 = PatchCollection(wedges, cmap=plt.cm.hsv, alpha=0.1, facecolor='grey')  # 修改颜色
    ax.add_collection(collection1)

    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.axis('equal')
    plt.axis('off')
    plt.legend(loc='best', frameon=True)
    plt.savefig(f'raresub_3method_circle_{subtype}.png', dpi=300, bbox_inches='tight')  # 增加dpi和bbox_inches参数，让图像更清晰
    plt.savefig(f'raresub_3method_circle_{subtype}.pdf', dpi=300, bbox_inches='tight')  # 增加dpi和bbox_inches参数，让图像更清晰
    plt.show()
