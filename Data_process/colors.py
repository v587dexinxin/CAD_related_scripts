# -*- coding: utf-8 -*-
"""
Created on Sat May 23 15:34:23 2026

@author: lenovo
"""


import matplotlib.pyplot as plt

# 常用科研作图颜色
colors = {
    "Blue": "#4C78A8",
    "Green": "#54A24B",
    "Red": "#E45756",
    "Purple": "#9370DB",
    "Light purple": "#B39DDB",
    "Lavender": "#D1C4E9",
    "Magenta": "#C51B7D",
    "Light magenta": "#D65F9E",
    "Very light magenta": "#E78AC3",
    "Orange": "#F28E2B",
    "Yellow": "#ECA400",
    "Cyan": "#72B7B2",
    "Grey": "#BDBDBD",
    "Light grey": "#D9D9D9",
    "Dark grey": "#666666",
    'light red1': "#F4A6A6",
    'light red2': "#F7B6B2",
    'light red3': "#FF9999",
    'light red4': "#EFA3A3",
    'light red5': "#F08080"
}

fig, ax = plt.subplots(figsize=(8, 5))

for i, (name, hex_code) in enumerate(colors.items()):
    ax.barh(i, 1, color=hex_code)
    ax.text(
        1.05, i,
        f"{name}   {hex_code}",
        va="center",
        fontsize=11
    )

ax.set_xlim(0, 2.8)
ax.set_ylim(-0.5, len(colors) - 0.5)
ax.axis("off")
ax.set_title("Common Python color codes", fontsize=14)

plt.tight_layout()

plt.savefig(
    "python_color_code_palette.pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "python_color_code_palette.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

import matplotlib.pyplot as plt

# 常用科研作图颜色
colors = {
    "Blue": "#4C78A8",
    "Green": "#54A24B",
    "Red": "#E45756",
    "Purple": "#9370DB",
    "Light purple": "#B39DDB",
    "Lavender": "#D1C4E9",
    "Magenta": "#C51B7D",
    "Light magenta": "#D65F9E",
    "Very light magenta": "#E78AC3",
    "Orange": "#F28E2B",
    "Yellow": "#ECA400",
    "Cyan": "#72B7B2",
    "Grey": "#BDBDBD",
    "Light grey": "#D9D9D9",
    "Dark grey": "#666666",
    'light red1': "#F4A6A6",
    'light red2': "#F7B6B2",
    'light red3': "#FF9999",
    'light red4': "#EFA3A3",
    'light red5': "#F08080",
    'light blue1': "#A6CEE3",
    'light blue2': "#9ECAE1",
    'light blue3': "#B3D7F2",
    'light blue4': "#ADD8E6",
    'light blue5': "#87CEFA",
    'light blue6': "#8ECAE6",
    'light blue7': "#BBDDF2",
    'light blue8': "#C6DBEF",
    'light blue9': "#6BAED6",
    'light blue10': "#4C78A8"
}

fig, ax = plt.subplots(figsize=(8, 5))

for i, (name, hex_code) in enumerate(colors.items()):
    ax.barh(i, 1, color=hex_code)
    ax.text(
        1.05, i,
        f"{name}   {hex_code}",
        va="center",
        fontsize=11
    )

ax.set_xlim(0, 2.8)
ax.set_ylim(-0.5, len(colors) - 0.5)
ax.axis("off")
ax.set_title("Common Python color codes", fontsize=14)

plt.tight_layout()

# plt.savefig(
#     "python_color_code_palette.pdf",
#     dpi=300,
#     bbox_inches="tight"
# )

# plt.savefig(
#     "python_color_code_palette.png",
#     dpi=300,
#     bbox_inches="tight"
# )

plt.show()