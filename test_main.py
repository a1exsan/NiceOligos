#import os
#from nicegui import app, ui

#os.makedirs('static_images', exist_ok=True)
#app.add_static_files('/img', 'static_images')

# Скопируйте файл в папку static_images предварительно
#ii = ui.interactive_image('/img/number_oligos_plot_1.png', cross='red',
#                          on_mouse=lambda e: e.sender.set_content(f'''
#                                  <circle cx="{e.image_x}" cy="{e.image_y}" r="5" fill="orange" />
#                              '''),
#                          )

#ui.run()

import os
from nicegui import events, ui, app

wells = {}
for x, num in zip(range(130, 1300, 100), '1 2 3 4 5 6 7 8 9 10 11 12'.split(' ')):
    for y, symb in zip(range(130, 900, 100), 'A B C D E F G H'.split(' ')):
        wells[(x, y)] = (symb, num)

print(wells)

def get_coord(x, y):
    out = {'well': ('', ''), 'key': (0, 0)}
    for key in wells.keys():
        r2 = (x - key[0])**2 + (y - key[1])**2
        if r2 <= 45**2:
            out['well'] = wells[key]
            out['key'] = key
            break
    return out

#get_coord(0,0)

def mouse_handler(e: events.MouseEventArguments):
    #image.content += f'<circle cx="{e.image_x}" cy="{e.image_y}" r="30" fill="none" stroke="red" stroke-width="4" />'
    #highlight.content = f'<circle cx="{e.image_x}" cy="{e.image_y}" r="28" fill="yellow" opacity="0.5" />'
    #print(e)

    if e.type == 'mousedown':
        coord = get_coord(e.image_x, e.image_y)
        if coord['well'] != ('', ''):
            highlight_1.content = f"<circle cx='{coord['key'][0]}' cy='{coord['key'][1]}' r='45' fill='orange' stroke='lightgreen' stroke-width='6' />"

    if e.type == 'mousemove':
        coord = get_coord(e.image_x, e.image_y)
        if coord['well'] != ('', ''):
            highlight_1.content = f"<circle cx='{coord['key'][0]}' cy='{coord['key'][1]}' r='45' fill='yellow' stroke='lightgreen' stroke-width='6' />"

os.makedirs('static_images', exist_ok=True)
app.add_static_files('/img', 'static_images')

image = ui.interactive_image('/img/plate_bkg_2.png', on_mouse=mouse_handler, #cross='red',
                             events=['mousedown', 'mousemove', 'mouseup', 'click'])

highlight_1 = image.add_layer()
highlight = image.add_layer()

for x in range(130, 1300, 100):
    for y in range(130, 900, 100):
        image.content += f'<circle cx="{x}" cy="{y}" r="45" fill="gray" stroke="lightgreen" stroke-width="6" />'

for x, num in zip(range(115, 1300, 100), '1 2 3 4 5 6 7 8 9 10 11 12'.split(' ')):
    image.content += f'<text x="{x}" y="{70}" fill="orange" font-size="48">{num}</text>'
    highlight.content += f'<text x="{x}" y="{920}" fill="lightblue" font-size="48">{num}</text>'

for y, num in zip(range(145, 900, 100), 'A B C D E F G H'.split(' ')):
    image.content += f'<text x="{40}" y="{y}" fill="orange" font-size="48">{num}</text>'
    highlight.content += f'<text x="{1300}" y="{y}" fill="lightblue" font-size="48">{num}</text>'


#ui.run()

