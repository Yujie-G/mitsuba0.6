## Biplane - Mitsuba renderer part

Forked from [Jiahui Fan's implement](https://github.com/sssssy/mitsuba-tensorNLB)


## Install on Ubuntu

### dependencies

```shell
sudo apt-get update
sudo apt-get install build-essential scons mercurial qt4-dev-tools libpng-dev libjpeg-dev libilmbase-dev libxerces-c-dev libboost-all-dev libopenexr-dev libglewmx-dev libxxf86vm-dev libpcrecpp0v5 libeigen3-dev libfftw3-dev
```



```shell
sudo add-apt-repository ppa:rock-core/qt4
sudo apt update
```

安装OpenEXR
  1. `pip install openexr`
  2. import error: https://zhuanlan.zhihu.com/p/615111375

### scons compile

build/Sconscript.install    line54-56 can be noted 

if QtWidgets was not found:

https://github.com/mitsuba-renderer/mitsuba/issues/125
 

```shell
scons -j 8
```

### add to bashrc

add:
```
export MITSUBA_PYVER=3.x
source /path/to/your/setpath.sh
```



## Known issues

-  [ubuntu20.04]: If **Unable to locate package qt4-dev-tools**:

    The Qt4 framework has been removed from Ubuntu 20.04 main repository.
    You can still get Qt4 libraries, adding the PPA rock-core/qt4
    Run in a terminal:

- 如果你的Python环境是3.8 && >=3.8.10，可能出现**libp11-kit.so.0: undefined symbol: ffi_type_pointer**,参考[这篇](https://blog.csdn.net/qq_38606680/article/details/129118491),降级Python到3.8.1

- 安装出的dist目录下没有Python插件.在`data\scons\detect_python.py`中line41-57中,找到相应目录,手动建立链接.`ln xxx.so sss.so`



