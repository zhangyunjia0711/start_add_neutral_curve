---
project: START_ALL_IN_ONE
summary: STability Analysis and Research Toolkit For LST, LPSE, NPSE and other methods.
project_github: none
project_download: none
author: Jianxin Liu@Tianjin University
author_pic: none
author_email: shookware@tju.edu.cn
src_dir: ../src
output_dir: ../doc/html
page_dir: ../doc/pages
exclude_dir: ../bin
             ../build
             ../src/third_party/cfgio
             ../src/include
exclude: functional.f90
         module_sort.f90
display: public
         protected
source: true
graph: true
search: true
proc_internals: true
print_creation_date: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            petsc:https://petsc.org
            fortran_stdlib:https://stdlib.fortran-lang.org
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
ordered_subpage: theory
ordered_subpage: guide
ordered_subpage: examples
css: ./css/style.css

---
[TOC]

介绍
========
**START_ALL_IN_ONE**是**START**套件中的一个组件。该程序主要是用来进行当地和全局的稳定性分析的。其中主要集成了线性稳定性理论分析(LST)，线性抛物化稳定性方程(LPSE)，非线性抛物化稳定性方程(NPSE)以及再次接触上沟通各个功能完成的各种进阶稳定性分析工具。

需求
========
代码安装需要支持Fortran2008以上标准的编译器，推荐Gfortran9.0+或intel fortran 2019+。同时代码还需要第三方库PETSc和SLEPc以及Fortran标准库fortran_stdlib。请在安装本代码之前先安装和配置如上软件。

文档
===
本代码文档采用Ford来自动生成。请预先安装ford，然后直接在项目目录下输入如下命令即可：
```bash
ford doc/docs.md
```

主要功能
========
本项目是一个开发中的项目，目前主要完成的功能和未来待开发的功能如下

- [X] LST(linear stability theory)
    - [X] Temporal mode
    - [X] spatial mode
    - [X] manual analysis with parameter variation
- [X] LPSE(linear parabolized stability equations)
    - [X] original/adjoint parabolized stability equations
- [X] NPSE(nonlinear parabolized stability equations)
- [ ] EPSE(expansion of the parabolized stability equation)
- [ ] Bi-global instability analysis in streamwise and spanwise directions
- [X] others
    - [X] optimal growth streak evaluation
    - [ ] parallel/nonparallel receptivity analysis
    - [ ] neutral curve generation
    - [ ] other features to be supplemented
- [X] friendly config file
- [X] parallel input/output (HDF5 binary files)
- [X] 文档化

License
===
本代码仅作为科研学习使用。或取代码请发邮件联系作者。使用本代码请引用作者相关工作。

商业应用授权相关事宜请联系作者shookware@tju.edu.cn详谈(刘建新@天津大学)