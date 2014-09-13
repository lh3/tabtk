Introduction
------------

Tabtk is a fast and lightweight tool for processing TAB-delimited formats.

Tabtk Examples
--------------

* Basic Unix `cut` (duplicated columns ignored):

        tabtk cut -f 5,1-3,6,6- file.txt

* Reorder columns:

        tabtk cut -rf 5,1-3,6 file.txt

* Duplicate columns (duplicated columns not ignored with option `-r`):

        tabtk cut -rf 1,1,1 file.txt

* Use both SPACE and TAB as the delimitor:

        tabtk cut -d space -f 1-3 file.txt

