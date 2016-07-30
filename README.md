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

* Cut a CSV file:

		tabtk cut -d csv -f 2-4 file.csv

  Commas can appear in double-quotation marks.

* Print lines in `streamed.txt` that matching `loaded.txt` on the first column:

		tabtk isct loaded.txt streamed.txt

* Print lines matching the first two columns:

		tabtk isct -1 1,2 loaded.txt streamed.txt

* Fixed-width view of a TAB delimited file and truncate long fields to 20

		tabtk view -l 20 tab.txt | less -S

  This by default loads 16MB data to RAM, not the whole file.

* Grep a pattern in specified columns:

		tabtk grep -f 2 "^rs[0-9]+" file.vcf

* Compute the mean, min and max of a numeric column:

		tabtk num -c 2 file.txt

* Compute the standard deviation and quartile:

		tabtk num -Qc2 file.txt
