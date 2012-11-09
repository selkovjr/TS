#!/bin/bash
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

run ()
{
    local EXIT_CODE=0
    eval $* >&2 || EXIT_CODE="$?"
    if [ ${EXIT_CODE} != 0 ]; then
        echo -e "ERROR: Status code '${EXIT_CODE}' while running\n\$ $*" >&2
        rm -f "${TSP_FILEPATH_PLUGIN_DIR}/${HTML_RESULTS}" "$JSON_RESULTS"
        exit 1
    fi
}
export -f run

write_html_header ()
{
    local HTML="${TSP_FILEPATH_PLUGIN_DIR}/header"
    if [ -n "$1" ]; then
        HTML="$1"
    fi
    local REFRESHRATE=0
    if [ -n "$2" ]; then
	REFRESHRATE=$2
    fi
    echo '<?xml version="1.0" encoding="iso-8859-1"?>' > "$HTML"
    echo '<!DOCTYPE html>' >> "$HTML"
    echo '<html>' >> "$HTML"
    print_html_head $REFRESHRATE >> "$HTML"
    echo '<title>Torrent Combine Barcodes Report</title>' >> "$HTML"
    echo '<style> body { background: #F0F0F4; padding: 0; margin: 8; } </style>' >> "$HTML"
    echo '<body>' >> "$HTML"
    print_html_logo >> "$HTML";

    if [ -z "$COV_PAGE_WIDTH" ];then
	echo '<div id="inner">' >> "$HTML"
    else
	echo "<div style=\"width:${COV_PAGE_WIDTH}px;margin-left:auto;margin-right:auto;height:100%\">" >> "$HTML"
    fi
    local OPTIONS=""
    if [ "$PLUGINCONFIG__COLLAPSE" = "Yes" ]; then
        OPTIONS="   Collapse Grouping"
    fi
    if [ "$PLUGINCONFIG__UNIQUESTARTS" = "Yes" ]; then
      OPTIONS="$OPTIONS   Create Unique Starts   Read End Resolution='$PLUGINCONFIG__LENGTHVAR'"
      if [ "$PLUGINCONFIG__GROUPUNIQUESTARTS" = "Yes" ]; then
          OPTIONS="$OPTIONS   Group Unique Starts"
      fi
    fi
    if [ "$PLUGINCONFIG__COMBINEGROUPS" = "Yes" ]; then
        OPTIONS="$OPTIONS   Combine Groups"
    fi
    echo "<h1><center><span style=\"cursor:help\" title=\"Barcode Grouping='${PLUGINCONFIG__GROUPING}'   Barcode Spacing='${PLUGINCONFIG__SPACING}'   Mapped Reads Threshold='${PLUGINCONFIG__READSTHRESHOLD}'${OPTIONS}\">Combine Barcodes Report</span></center></h1>" >> "$HTML"
}

write_html_footer ()
{
    local HTML="${TSP_FILEPATH_PLUGIN_DIR}/footer"
    if [ -n "$1" ]; then
        HTML="$1"
    fi
    print_html_end_javascript >> "$HTML"
    print_html_footer >> "$HTML"
    echo ' <br/><br/></div>' >> "$HTML"
    echo '</body></html>' >> "$HTML"
}

display_static_progress ()
{
    local HTML="${TSP_FILEPATH_PLUGIN_DIR}/${HTML_RESULTS}"
    if [ -n "$1" ]; then
        HTML="$1"
    fi
    echo "<br/><h3 style=\"text-align:center;color:red\">*** Processing is not complete ***</h3>" >> "$HTML"
    echo "<a href=\"javascript:document.location.reload();\" ONMOUSEOVER=\"window.status='Refresh'; return true\">" >> "$HTML"
    echo "<div style=\"text-align:center\">Click here to refresh</div></a>" >> "$HTML"
}

write_json ()
{
    echo "{" > "$JSON_RESULTS"
    echo "  \"Group Name\" : \"$PLUGIN_RESULTS_NAME\"," >> "$JSON_RESULTS"
    echo "  \"Group Size\"  : \"$PLUGIN_NUM_RUNS\"," >> "$JSON_RESULTS"
    echo "  \"Members\"  : \"$PLUGIN_RUN_LIST\"," >> "$JSON_RESULTS"
    echo "  \"Total AQ17\" : \"$PLUGIN_TOTAL_AQ17\"" >> "$JSON_RESULTS"
    echo "}" >> "$JSON_RESULTS"
}

