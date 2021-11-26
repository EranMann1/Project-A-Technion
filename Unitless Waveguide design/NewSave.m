if exist('Title','var')
else
  Title = 'random title';
end
string = [pwd '\Saved Params\' date '\' Title '.mat'];
if exist('Saved Params','dir')
else
    mkdir('Saved Params');
end
cd([pwd '\Saved Params']);
if exist(date,'dir')
else
    mkdir(date);
end
cd([pwd '\' date]);
if exist(Title,'file') 
    warndlg('title already exsists today');
else
  savefig(a1,[Title '_1.fig']); 
  savefig(a2,[Title '_2.fig']);
  savefig(a3,[Title '_3.fig']);
  savefig(a4,[Title '_4.fig']);
  clear a1 a2 a3 a4;
  save(Title);
end
cd('..');
cd('..');