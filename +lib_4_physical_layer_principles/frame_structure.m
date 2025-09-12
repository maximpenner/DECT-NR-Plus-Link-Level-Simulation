function [T_frame, N_FRAME_slot, T_slot] = frame_structure()
    T_frame = 10e-3;
    N_FRAME_slot = 24;
    T_slot = T_frame/N_FRAME_slot;
end
