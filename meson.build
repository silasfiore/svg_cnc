project('svg_cnc','c')

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
x11_dep = dependency('X11')

deps = [m_dep,x11_dep]

build_target = ['src/main.c',
                'src/svg_path.c',
                'src/graphic.c',
                'src/bezdist.c',
                'src/qbez.c']
include_dirs = 'inc'

executable('prog',build_target,include_directories : include_dirs,dependencies : deps)
executable('dat2svg','src/dat2svg.c',dependencies : [m_dep])