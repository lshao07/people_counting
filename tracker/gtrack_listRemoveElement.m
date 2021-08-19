function list = gtrack_listRemoveElement(list, elem)
    index = find(list.list == elem); 
    list.list(index) = [];
    list.count = list.count - 1;
end