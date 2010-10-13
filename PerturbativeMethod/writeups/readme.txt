# Start: Sept 2010

Detailed notes for filling in steps in the derivations of the Alard
2007 and 2008 papers, so far.


Main File: Report_on_Perturbative_Method.tex


Basically the content of the file for now:

-title.tex : This file has the title, author and institution.

-introduction.tex: Must be contain the introduction of the reports, if someone can edit, thanks.

-perturbative_method_chap2.tex (notes_forAlard2007Paper.tex)

-perturbative_method_chap3.tex (notes_forAlard2008Paper.tex)

************************************************************************************************
* If you want add section or subsection, only edit these files.

* If you want add some chapter, create the file

perturbative_method_chap#.tex (where # is the number of chapter) 
and use the command 
\include{perturbative_method_chap#} (inside the Report_on_Perturbative_Method.tex file)

******************************************************************************************************
Format graphics accepted: *.jpg, *.png, *.pdf (doesn't accept *ps or *eps)
The file graphics must be in the graphics folder, and for include them use the typical latex figure commands.
If you want edit some figure, use the GIMP, and save the final file in *.png (for remain the resolution)


1.- For Wrap Figure (figure inside the text)

\begin{wrapfigure}{r}{0.4\textwidth}
  \begin{center}
   \includegraphics[width=0.35\textwidth]{graphics/sourceplane.pdf}
  \end{center}
    \caption{\label{circular_source}Circular source in the source plane.}
\end{wrapfigure}

Note 1) {r} indicates the position of the figure, right or left side
Note 2) The width in includegraphics must be less than \textwidth.

 
2.- For single figure.

\begin{figure}[position]
\includegraphics[dimensions of the figure]{graphics/file_name.format}
\caption{\label{} Caption here}
\end{figure}


3.- For two figures

\begin{figure*}[!htp]
\resizebox{\hsize}{!}{
\subfigure{\includegraphics{graphics/file_left_name.format}}
\hspace{3.cm} ! you can change this value...
\subfigure{\includegraphics{graphics/file_right.format}}}
\caption{\label{measure_arcs} Caption here}
\end{figure*}

* I recomend see the subfigure package manual, to include more than 2 figures.
 * The comand \hsize, adjust the size of figures to the page width, some doubts, read the A&A manual.

******************************************************************************************************************88
See the user definitions, for example

\def \re {R_{_\mathrm{E}}}

You are free to include new definitions.

Please, use the typical command for fraction \frac{expression1}{expression2} or \dfrac{expression1}{expression2}, 
not: expression_1 \over expression_2

********************************************************************************************************************
To compile: pdflatex Report_on_Perturbative_Method.tex