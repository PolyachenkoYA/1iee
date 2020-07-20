#!/usr/bin/gawk -f
function add_res(list, line){
    match(line, /#([0-9]+)\:/, arr)
    if(arr[1] < 8){
        match(line, /(#[0-9]+\:[0-9]+)\./, arr)
        res_id = arr[1]
        return ((res_id in list) ? "0" : res_id)
    }
}
function print_res(line, printed){
    match(line, /#([0-9]+)\:/, arr)
    if(arr[1] < 8){
        match(line, /(#[0-9]+\:[0-9]+)\./, arr)
        if(printed == 1){
            printf ","
        }
        printf "%s",arr[1]
        printed = 1
    }
    return printed
}
BEGIN{
    printed = 0;
}
{
    if(($4 < 0.01) && ($4 != "")){
        printed = print_res($1, printed)
        printed = print_res($2, printed)
    }
}
END{
    print " "
}
