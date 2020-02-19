function cost_out = costgen1(s)
    cost_out =0; 

    S1 = stepinfo(s,'SettlingTimeThreshold',0.01);
    cost_out = S1.SettlingTime +S1.Overshoot/100+S1.Undershoot/100 +abs(1-S1.Peak);

end

