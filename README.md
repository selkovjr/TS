# TVC candidate generator

This fork of TS was made to isolate candidate generator functionality (originally borrowed from `freebayes`) from the rest of Ion Torrent stuff. It is the fastest way to call raw variants in a BAM.

It depends on Ion Torrent's version of `bamtools`, which it pulls from [http://updates.iontorrent.com/updates/software/external/bamtools-2.4.0.20150702+git15eadb925f.tar.gz](http://updates.iontorrent.com/updates/software/external/bamtools-2.4.0.20150702+git15eadb925f.tar.gz)

It also pulls specific verions of `htslib`, `samtools`, and `picard-tools` from the same site.

Original build instructions:

> Please see buildTools/BUILD.txt for build requirements and instructions.
> Note especially the section BUILD SPECIFIC MODULES.
>
> To build Analysis module, the command is:
> MODULES=Analysis ./buildTools/build.sh

Because this repo has been trimmed to a small portion of TVC's original code, which is part of the Analysis module, the build script has been updated to build that by default:

```
rm -rf build;  ./buildTools/build.sh
```
