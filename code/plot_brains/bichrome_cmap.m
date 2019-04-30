function cmap = bichrome_cmap(min_col,max_col)
    length = 100;    
    cmap = [linspace(min_col(1),max_col(1),length)', linspace(min_col(2),max_col(2),length)', linspace(min_col(3),max_col(3),length)'];