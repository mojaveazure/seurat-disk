% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loom.R
\docType{class}
\name{loom-class}
\alias{loom-class}
\alias{loom}
\title{A class for connections to loom files}
\format{
An \code{\link[R6]{R6Class}} object
}
\description{
A class for connections to loom files

A class for connections to loom files
}
\seealso{
\code{\link[hdf5r]{H5File}}
}
\section{Super classes}{
\code{\link[hdf5r:H5RefClass]{hdf5r::H5RefClass}} -> \code{\link[hdf5r:H5File]{hdf5r::H5File}} -> \code{\link[SeuratDisk:scdisk]{SeuratDisk::scdisk}} -> \code{loom}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-loom-add_attribute}{\code{loom$add_attribute()}}
\item \href{#method-loom-add_graph}{\code{loom$add_graph()}}
\item \href{#method-loom-add_layer}{\code{loom$add_layer()}}
\item \href{#method-loom-version}{\code{loom$version()}}
\item \href{#method-loom-timestamp}{\code{loom$timestamp()}}
\item \href{#method-loom-last.modified}{\code{loom$last.modified()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5RefClass" data-id="close"><a href='../../hdf5r/html/H5RefClass.html#method-H5RefClass-close'><code>hdf5r::H5RefClass$close()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5RefClass" data-id="dec_ref"><a href='../../hdf5r/html/H5RefClass.html#method-H5RefClass-dec_ref'><code>hdf5r::H5RefClass$dec_ref()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5RefClass" data-id="get_file_id"><a href='../../hdf5r/html/H5RefClass.html#method-H5RefClass-get_file_id'><code>hdf5r::H5RefClass$get_file_id()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5RefClass" data-id="get_obj_type"><a href='../../hdf5r/html/H5RefClass.html#method-H5RefClass-get_obj_type'><code>hdf5r::H5RefClass$get_obj_type()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5RefClass" data-id="get_ref"><a href='../../hdf5r/html/H5RefClass.html#method-H5RefClass-get_ref'><code>hdf5r::H5RefClass$get_ref()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5RefClass" data-id="inc_ref"><a href='../../hdf5r/html/H5RefClass.html#method-H5RefClass-inc_ref'><code>hdf5r::H5RefClass$inc_ref()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5RefClass" data-id="methods"><a href='../../hdf5r/html/H5RefClass.html#method-H5RefClass-methods'><code>hdf5r::H5RefClass$methods()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_delete"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_delete'><code>hdf5r::H5File$attr_delete()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_delete_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_delete_by_idx'><code>hdf5r::H5File$attr_delete_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_delete_by_name"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_delete_by_name'><code>hdf5r::H5File$attr_delete_by_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_exists"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_exists'><code>hdf5r::H5File$attr_exists()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_exists_by_name"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_exists_by_name'><code>hdf5r::H5File$attr_exists_by_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_get_number"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_get_number'><code>hdf5r::H5File$attr_get_number()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_info_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_info_by_idx'><code>hdf5r::H5File$attr_info_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_info_by_name"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_info_by_name'><code>hdf5r::H5File$attr_info_by_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_name_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_name_by_idx'><code>hdf5r::H5File$attr_name_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_open"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_open'><code>hdf5r::H5File$attr_open()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_open_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_open_by_idx'><code>hdf5r::H5File$attr_open_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_open_by_name"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_open_by_name'><code>hdf5r::H5File$attr_open_by_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_rename"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_rename'><code>hdf5r::H5File$attr_rename()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="attr_rename_by_name"><a href='../../hdf5r/html/H5File.html#method-H5File-attr_rename_by_name'><code>hdf5r::H5File$attr_rename_by_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="close_all"><a href='../../hdf5r/html/H5File.html#method-H5File-close_all'><code>hdf5r::H5File$close_all()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="commit"><a href='../../hdf5r/html/H5File.html#method-H5File-commit'><code>hdf5r::H5File$commit()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="create_attr"><a href='../../hdf5r/html/H5File.html#method-H5File-create_attr'><code>hdf5r::H5File$create_attr()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="create_attr_by_name"><a href='../../hdf5r/html/H5File.html#method-H5File-create_attr_by_name'><code>hdf5r::H5File$create_attr_by_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="create_dataset"><a href='../../hdf5r/html/H5File.html#method-H5File-create_dataset'><code>hdf5r::H5File$create_dataset()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="create_group"><a href='../../hdf5r/html/H5File.html#method-H5File-create_group'><code>hdf5r::H5File$create_group()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="create_reference"><a href='../../hdf5r/html/H5File.html#method-H5File-create_reference'><code>hdf5r::H5File$create_reference()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="exists"><a href='../../hdf5r/html/H5File.html#method-H5File-exists'><code>hdf5r::H5File$exists()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="file_info"><a href='../../hdf5r/html/H5File.html#method-H5File-file_info'><code>hdf5r::H5File$file_info()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="flush"><a href='../../hdf5r/html/H5File.html#method-H5File-flush'><code>hdf5r::H5File$flush()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="get_filename"><a href='../../hdf5r/html/H5File.html#method-H5File-get_filename'><code>hdf5r::H5File$get_filename()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="get_filesize"><a href='../../hdf5r/html/H5File.html#method-H5File-get_filesize'><code>hdf5r::H5File$get_filesize()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="get_intent"><a href='../../hdf5r/html/H5File.html#method-H5File-get_intent'><code>hdf5r::H5File$get_intent()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="get_obj_count"><a href='../../hdf5r/html/H5File.html#method-H5File-get_obj_count'><code>hdf5r::H5File$get_obj_count()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="get_obj_ids"><a href='../../hdf5r/html/H5File.html#method-H5File-get_obj_ids'><code>hdf5r::H5File$get_obj_ids()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="get_obj_name"><a href='../../hdf5r/html/H5File.html#method-H5File-get_obj_name'><code>hdf5r::H5File$get_obj_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="group_info"><a href='../../hdf5r/html/H5File.html#method-H5File-group_info'><code>hdf5r::H5File$group_info()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="group_info_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-group_info_by_idx'><code>hdf5r::H5File$group_info_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="group_info_by_name"><a href='../../hdf5r/html/H5File.html#method-H5File-group_info_by_name'><code>hdf5r::H5File$group_info_by_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link"><a href='../../hdf5r/html/H5File.html#method-H5File-link'><code>hdf5r::H5File$link()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_copy_from"><a href='../../hdf5r/html/H5File.html#method-H5File-link_copy_from'><code>hdf5r::H5File$link_copy_from()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_copy_to"><a href='../../hdf5r/html/H5File.html#method-H5File-link_copy_to'><code>hdf5r::H5File$link_copy_to()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_create_external"><a href='../../hdf5r/html/H5File.html#method-H5File-link_create_external'><code>hdf5r::H5File$link_create_external()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_create_hard"><a href='../../hdf5r/html/H5File.html#method-H5File-link_create_hard'><code>hdf5r::H5File$link_create_hard()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_create_soft"><a href='../../hdf5r/html/H5File.html#method-H5File-link_create_soft'><code>hdf5r::H5File$link_create_soft()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_delete"><a href='../../hdf5r/html/H5File.html#method-H5File-link_delete'><code>hdf5r::H5File$link_delete()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_delete_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-link_delete_by_idx'><code>hdf5r::H5File$link_delete_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_exists"><a href='../../hdf5r/html/H5File.html#method-H5File-link_exists'><code>hdf5r::H5File$link_exists()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_info"><a href='../../hdf5r/html/H5File.html#method-H5File-link_info'><code>hdf5r::H5File$link_info()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_info_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-link_info_by_idx'><code>hdf5r::H5File$link_info_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_move_from"><a href='../../hdf5r/html/H5File.html#method-H5File-link_move_from'><code>hdf5r::H5File$link_move_from()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_move_to"><a href='../../hdf5r/html/H5File.html#method-H5File-link_move_to'><code>hdf5r::H5File$link_move_to()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_name_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-link_name_by_idx'><code>hdf5r::H5File$link_name_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_value"><a href='../../hdf5r/html/H5File.html#method-H5File-link_value'><code>hdf5r::H5File$link_value()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="link_value_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-link_value_by_idx'><code>hdf5r::H5File$link_value_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="ls"><a href='../../hdf5r/html/H5File.html#method-H5File-ls'><code>hdf5r::H5File$ls()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="mount"><a href='../../hdf5r/html/H5File.html#method-H5File-mount'><code>hdf5r::H5File$mount()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="obj_copy_from"><a href='../../hdf5r/html/H5File.html#method-H5File-obj_copy_from'><code>hdf5r::H5File$obj_copy_from()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="obj_copy_to"><a href='../../hdf5r/html/H5File.html#method-H5File-obj_copy_to'><code>hdf5r::H5File$obj_copy_to()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="obj_info"><a href='../../hdf5r/html/H5File.html#method-H5File-obj_info'><code>hdf5r::H5File$obj_info()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="obj_info_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-obj_info_by_idx'><code>hdf5r::H5File$obj_info_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="obj_info_by_name"><a href='../../hdf5r/html/H5File.html#method-H5File-obj_info_by_name'><code>hdf5r::H5File$obj_info_by_name()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="open"><a href='../../hdf5r/html/H5File.html#method-H5File-open'><code>hdf5r::H5File$open()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="open_by_idx"><a href='../../hdf5r/html/H5File.html#method-H5File-open_by_idx'><code>hdf5r::H5File$open_by_idx()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="path_valid"><a href='../../hdf5r/html/H5File.html#method-H5File-path_valid'><code>hdf5r::H5File$path_valid()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="print"><a href='../../hdf5r/html/H5File.html#method-H5File-print'><code>hdf5r::H5File$print()</code></a></span></li>
<li><span class="pkg-link" data-pkg="hdf5r" data-topic="H5File" data-id="unmount"><a href='../../hdf5r/html/H5File.html#method-H5File-unmount'><code>hdf5r::H5File$unmount()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeuratDisk" data-topic="scdisk" data-id="chunk.points"><a href='../../SeuratDisk/html/scdisk.html#method-scdisk-chunk.points'><code>SeuratDisk::scdisk$chunk.points()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeuratDisk" data-topic="scdisk" data-id="finalizer"><a href='../../SeuratDisk/html/scdisk.html#method-scdisk-finalizer'><code>SeuratDisk::scdisk$finalizer()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SeuratDisk" data-topic="scdisk" data-id="initialize"><a href='../../SeuratDisk/html/scdisk.html#method-scdisk-initialize'><code>SeuratDisk::scdisk$initialize()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-loom-add_attribute"></a>}}
\if{latex}{\out{\hypertarget{method-loom-add_attribute}{}}}
\subsection{Method \code{add_attribute()}}{
Add an attribute
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{loom$add_attribute(x, name, type = c("global", "row", "col"))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{Object to add as an attribute}

\item{\code{name}}{Name to store attribute as}

\item{\code{type}}{Type of attribute to add}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-loom-add_graph"></a>}}
\if{latex}{\out{\hypertarget{method-loom-add_graph}{}}}
\subsection{Method \code{add_graph()}}{
Add a graph
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{loom$add_graph(x, name, type = c("col", "row"), verbose = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{...}

\item{\code{name}}{...}

\item{\code{type}}{...}

\item{\code{verbose}}{...}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-loom-add_layer"></a>}}
\if{latex}{\out{\hypertarget{method-loom-add_layer}{}}}
\subsection{Method \code{add_layer()}}{
Add a layer to this loom file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{loom$add_layer(x, name, transpose = TRUE, verbose = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{An object to save as a layer}

\item{\code{name}}{Name to store layer as}

\item{\code{transpose}}{...}

\item{\code{verbose}}{...}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Invisibly returns \code{NULL}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-loom-version"></a>}}
\if{latex}{\out{\hypertarget{method-loom-version}{}}}
\subsection{Method \code{version()}}{
Get version information
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{loom$version()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link[base]{numeric_version}} object with the loom
specification version information
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-loom-timestamp"></a>}}
\if{latex}{\out{\hypertarget{method-loom-timestamp}{}}}
\subsection{Method \code{timestamp()}}{
Add a timestamp to a dataset or group as an HDF5 attribute
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{loom$timestamp(name = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{Name of dataset or group to add timestamp to; if \code{NULL},
timestamps the file as a whole}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Invisibly returns the object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-loom-last.modified"></a>}}
\if{latex}{\out{\hypertarget{method-loom-last.modified}{}}}
\subsection{Method \code{last.modified()}}{
Retrieve a timestamp from a dataset or group
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{loom$last.modified(name = NULL, locale = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{Name of dataset or group to retrieve timestamp from; if
\code{NULL}, retrieves timestamp from at the file-level}

\item{\code{locale}}{Change the timestamp of to the timezone of the locale}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character with the timestamp
}
}
}
