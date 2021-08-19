function [elem, list] = gtrack_listDequeue(list)
    list.count = list.count - 1;
    elem = list.list(1);
    list.list(1) = [];
end