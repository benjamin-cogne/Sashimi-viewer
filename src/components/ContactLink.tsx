/**
 * A "Contact" link that works without a mail program: a mailto: link does nothing where none is set up (a webmail
 * only, a browser without a mail handler). It opens a small card instead: the address, to read or copy, and a new
 * message in the mail program, in Outlook on the web or in Gmail.
 */
import { useEffect, useRef, useState } from 'react';

export default function ContactLink({ email, subject, label = 'Contact', up = false, className = '' }: {
  email: string; subject: string; label?: string;
  /** open the card above the link (a footer at the bottom of the page) */
  up?: boolean;
  className?: string;
}) {
  const [open, setOpen] = useState(false);
  const [copied, setCopied] = useState(false);
  const box = useRef<HTMLSpanElement>(null);
  useEffect(() => {
    if (!open) return;
    const away = (e: MouseEvent) => { if (box.current && !box.current.contains(e.target as Node)) setOpen(false); };
    const esc = (e: KeyboardEvent) => { if (e.key === 'Escape') { e.stopPropagation(); setOpen(false); } };
    document.addEventListener('mousedown', away);
    window.addEventListener('keydown', esc, true);
    return () => { document.removeEventListener('mousedown', away); window.removeEventListener('keydown', esc, true); };
  }, [open]);
  const to = encodeURIComponent(email), su = encodeURIComponent(subject);
  const copy = () => {
    const done = () => { setCopied(true); setTimeout(() => setCopied(false), 1500); };
    navigator.clipboard?.writeText(email).then(done, () => { /* blocked: the address stays selectable */ });
  };
  const item = 'block px-2.5 py-1 rounded text-[12px] text-slate-700 hover:bg-slate-100 no-underline';
  return (
    <span ref={box} className="relative inline-block">
      <button type="button" onClick={() => setOpen(o => !o)} aria-expanded={open} title={`Write to ${email}`}
        className={`underline decoration-dotted ${className}`}>{label}</button>
      {open && (
        <span role="dialog" aria-label={`Contact ${email}`}
          className={`absolute z-[90] ${up ? 'bottom-full mb-1.5' : 'top-full mt-1.5'} right-0 w-[260px] rounded-lg border border-slate-200 bg-white shadow-lg p-2 text-left text-slate-800 normal-case tracking-normal`}>
          <span className="flex items-center gap-1.5 px-1">
            <span className="flex-1 font-mono text-[12px] select-all break-all">{email}</span>
            <button type="button" onClick={copy} className="shrink-0 text-[11px] px-1.5 py-0.5 rounded border border-slate-300 text-slate-600 hover:bg-slate-50">{copied ? '✓ copied' : 'Copy'}</button>
          </span>
          <span className="block mt-1.5 border-t border-slate-100 pt-1">
            <a className={item} href={`mailto:${email}?subject=${su}`} onClick={() => setOpen(false)}>Open in my mail program</a>
            <a className={item} href={`https://outlook.office.com/mail/deeplink/compose?to=${to}&subject=${su}`} target="_blank" rel="noopener noreferrer" onClick={() => setOpen(false)}>Write in Outlook (web)</a>
            <a className={item} href={`https://mail.google.com/mail/?view=cm&fs=1&to=${to}&su=${su}`} target="_blank" rel="noopener noreferrer" onClick={() => setOpen(false)}>Write in Gmail</a>
          </span>
          <span className="block px-2.5 pt-1 text-[10.5px] text-slate-400">If nothing opens, copy the address.</span>
        </span>
      )}
    </span>
  );
}
