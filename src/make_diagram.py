import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots(1, 1, figsize=(6, 6))

square = patches.Rectangle((0, 0), 1, 1, linewidth=3, edgecolor='black', facecolor='white')
ax.add_patch(square)

samples = [
    (0.3, 0.3, 0.25, 'bottom-left'),
    (0.35, 0.7, 0.40, 'top-left'),
    (0.6, 0.6, 0.25, 'top-right'),
    (0.7, 0.4, 0.10, 'bottom-right')
]

colors = ['#e74c3c', '#e74c3c', '#e74c3c', '#e74c3c']  
for i, (x, y, weight, label) in enumerate(samples):
    size = weight * 2000
    ax.scatter(x, y, s=size, c=[colors[i]], alpha=0.7, edgecolors='black', linewidth=2, zorder=3)
    
    ax.text(x, y - 0.08, f'({x}, {y}, {weight})', ha='center', va='top', fontsize=10, fontweight='bold')

ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-0.1, 1.1)
ax.set_aspect('equal')
ax.set_xlabel('x', fontsize=14, fontweight='bold')
ax.set_ylabel('y', fontsize=14, fontweight='bold')
ax.set_title('Custom Sample Pattern\n(4 samples with non-uniform weights)', fontsize=14, fontweight='bold', pad=20)

ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xticks([0, 0.5, 1])
ax.set_yticks([0, 0.5, 1])

ax.plot(0.5, 0.5, 'k+', markersize=15, markeredgewidth=2, label='Pixel Center (0.5, 0.5)')
ax.legend(loc='upper right', framealpha=0.9)

plt.tight_layout()
plt.savefig('/Users/abhinavmahajan/Documents/Sem2/Computer Graphics/Scotty3D/assignments/A1-writeup/student/task7-custom-pattern.png', dpi=150, bbox_inches='tight')
plt.show()

print("Sample pattern visualization saved as 'task7-custom-pattern.png'")