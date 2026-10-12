/**
 * Moving a floating panel (a site's panel, the phasing panel): pressed anywhere on it but its controls (fields, buttons,
 * menus, links) and its scroll bar, then moved. The panel follows once the pointer has gone a few pixels, so a plain
 * click still clicks; the moves are listened to on the window, so a fast drag that leaves the panel keeps it.
 */
import type React from 'react';

export interface PanelPlace { left: number; top: number }

/** Pixels the pointer moves before the panel follows. */
const SLOP = 4;

export function startPanelDrag(e: React.PointerEvent<HTMLElement>, at: PanelPlace, move: (p: PanelPlace) => void): void {
  if (e.button !== 0 || !e.isPrimary) return;
  const target = e.target as Element;
  if (target.closest('input, textarea, select, button, a, label, [data-no-drag]')) return;
  const el = e.currentTarget, box = el.getBoundingClientRect();
  // the scroll bar of a panel taller than the window
  if (el.scrollHeight > el.clientHeight && e.clientX > box.left + el.clientLeft + el.clientWidth) return;
  const x0 = e.clientX, y0 = e.clientY, width = box.width;
  let moving = false;
  const onMove = (ev: PointerEvent) => {
    if (ev.pointerId !== e.pointerId) return;
    if (!moving) {
      if (Math.hypot(ev.clientX - x0, ev.clientY - y0) < SLOP) return;
      moving = true;
      document.body.style.userSelect = 'none';
      window.getSelection()?.removeAllRanges();
    }
    ev.preventDefault();
    move({
      // kept within reach: at least 60 px of it on screen, its head never above the window
      left: Math.max(60 - width, Math.min(at.left + ev.clientX - x0, window.innerWidth - 60)),
      top: Math.max(0, Math.min(at.top + ev.clientY - y0, window.innerHeight - 40)),
    });
  };
  const onUp = (ev: PointerEvent) => {
    if (ev.pointerId !== e.pointerId) return;
    window.removeEventListener('pointermove', onMove);
    window.removeEventListener('pointerup', onUp);
    window.removeEventListener('pointercancel', onUp);
    document.body.style.userSelect = '';
  };
  window.addEventListener('pointermove', onMove);
  window.addEventListener('pointerup', onUp);
  window.addEventListener('pointercancel', onUp);
}
