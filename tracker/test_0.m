clear all;

array.count = 5;
array.list = [1, 2, 3, 4, 5];

element = 10;

[element_1,array] = gtrack_listDequeue(array);


array = gtrack_listEnqueue(array, element);

flag = gtrack_isListEmpty(array);

array = gtrack_listRemoveElement(array, 3);