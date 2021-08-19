function list = gtrack_listEnqueue(list, elem)
    list.count = list.count + 1;
    list.list(list.count) = elem;

end